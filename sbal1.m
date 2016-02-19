function [U,P,ress] = sbal1(U,P,u,iter, verbose)
% Sparse Bundle adjustment for calibrated cameras using the L1 norm
%
% inputs  - u: 1xD cell with image data.
%           u{i} is of size 3xN, where N is the number of observed points.
%           If point j is not observed in image then u{i}(:,j) = NaN.
%
%         - U: 3xD cell with 3D points
%
%         - P: 1xD cell with camera matrices
%
%         - iter: number of iterations
%
% outputs - U: 3xD cell with 3D points
%
%         - P: 1xD cell with camera matrices
%
% (C)  2011 Zhijun Dai (ISCAS) zhijun@iscas.ac.cn, dai.zhijun@gmail.com
%  The code of setup lin system is based on 2010 Carl Olsson (calle@maths.lth.se, carl.a.c.olsson@gmail.com)
%
if nargin <4
    iter = 50;
end
if nargin <5
    verbose = 1;
end

mu=1;
i=1;
stop=0;
maxiter=iter;
if size(U,1) == 3
    U = [U; ones(1,size(U,2))];
end
if verbose
    fprintf('Iter:\t Error L2:\t Error L1:\n');
end
res1 = compute_res_l1(P,U,u);
res2 = compute_res_l2(P,U,u);
if verbose
    fprintf('%d\t%f\t%f\n',0,res2,res1);
end
ress = [res1; res2];
while ~stop && i<=maxiter    
    [A,B, uStep,UrowIndex] = setup_lin_system(P,U,u);
    res1 = compute_res_l1(P,U,u);   
    % Solve first order approximation in U,
    %d = compute_update(A,-B);
    %d = compute_update_trust(A,-B,mu);
    d = compute_update_interior(A,-B);
    %d = compute_update_mosek(A,-B);
    %d = compute_update_sparse_l1(A,-B,UrowIndex);
    %d = compute_update_interior_trust(A,-B);
    %d = compute_update_trust2(A,-B, mu);
    %d = compute_update_simplex(A,-B);
    [Pnew,Unew] = update_var(d,P,U,u);
    res1new = compute_res_l1(Pnew,Unew,u);
    res2 = compute_res_l2(Pnew,Unew,u);
    
    gain= res1-res1new;
    while gain < 0 && norm(d,1) > 1e-7
        d = d/2;
        mu = mu/2;
        [Pnew,Unew] = update_var(d,P,U,u);
        res1new = compute_res_l1(Pnew,Unew,u);
        res2 = compute_res_l2(Pnew,Unew,u);
        gain= res1-res1new;
    end
    if gain > 0
        U = Unew;
        P = Pnew;
    end
    % --- Stopping criterion ---
    if (res1 - res1new) < 1e-6 || norm(d,1) < 1e-7
        stop=1;
    end
    if verbose
        fprintf('%d\t%f\t%f\n',i,res2,res1new);
    end
    ress=[ress [res1new; res2]];
    i=i+1;
end
U = U(1:3,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%algorithm use simplex method
function d = compute_update(A_in,B)
[m,n]=size(A_in);
f=[sparse(n,1);ones(m,1)];
%f=[zeros(m,1);ones(m,1)];
A=[A_in -speye(m,m);-A_in -speye(m,m)];
b=[B;-B];
x = linprog(f,A,b);
d = x(1:n);
function d = compute_update_mosek(A_in,B)
[m,n]=size(A_in);
f=[sparse(n,1);ones(m,1)];
%f=[zeros(m,1);ones(m,1)];
A=[A_in -speye(m,m);-A_in -speye(m,m)];
b=[B;-B];
x = linprog2(f,A,b);
d = x(1:n);

function d = compute_update_trust(J,B,mu)
opt=optimset; opt.Display='off';
[m,n]=size(J);
%f=[sparse(n,1);ones(m,1)];
%f=[zeros(m,1);ones(m,1)];
%A=[A_in -speye(m,m);-A_in -speye(m,m)];
%b=[B;-B];
f=[sparse(n,1); ones(m,1); sparse(n,1)];
A=[ [-J  -speye(m) sparse(m,n)];
    [ J  -speye(m) sparse(m,n)];
    [ speye(n) sparse(n,m) -speye(n)];
    [ -speye(n) sparse(n,m) -speye(n)];
    [ sparse(1,n) sparse(1,m) ones(1,n)]
    ];
b=[-B ; B ; zeros(n,1); zeros(n,1); mu];
[d,fval,exitflag,output,lambdaV]=linprog(f,A,b,[],[],[],[],[],opt);
d=d(1:n);

function d = compute_update_trust2(J,B,mu)
opt=optimset; opt.Display='off';
[m,n]=size(J);
%f=[sparse(n,1);ones(m,1)];
%f=[zeros(m,1);ones(m,1)];
%A=[A_in -speye(m,m);-A_in -speye(m,m)];
%b=[B;-B];
f=[sparse(n,1); ones(m,1)];
A=[ [-J  -speye(m)];
    [ J  -speye(m)];
    [ speye(n) sparse(n,m)];
    [ -speye(n) sparse(n,m)]
    ];
b=[-B ; B ; ones(n,1)*mu; ones(n,1)*mu];
[d,fval,exitflag,output,lambdaV]=linprog(f,A,b,[],[],[],[],[],opt);
d=d(1:n);

%function [V,z,basic]=compute_update_simplex(A,B)
function z=compute_update_simplex(A,B)
%
[m,n]=size(A);
% Form LP
fv=[sparse(n,1); ones(m,1); sparse(2*m,1)];
Av=[ [[-A -speye(m) ]; [A  -speye(m) ]] speye(2*m) ];
bv=[-B; B ];

% Construct starting point for simplex phase 2
v=zeros(n,1);
res=Av(1:m,1:n)*v-bv(1:m);
t=abs(res);
s=[t-res ; t+res];
x0=[v;t;s];
small=1e-6;
basic0=[find((t)>=small )]+n;
basic0=[basic0; [find((s)>=2*small )]+m+n];
ind=find(t<small);
basic0=[basic0 ; ind+n+m ; ind+n+m*2];
basic0=sort(basic0);
%[z2,basic2]=simplexfas2b(fv,Av,bv,m,r,n,W,[]);
[z,fval,basic]=simplexphase2(fv,Av,bv,x0,basic0);

% Output
basic=sort(basic);
%V=reshape(z(1:n*r)-z(n*r+1:2*n*r),[r n]);
z=z(1:n);
1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = compute_res_l2(P,U,u)
res = 0;
for i = 1:length(P);
    %vis = isfinite(u{i}(1,:));
    vis = find(u{i}(1,:)<inf);
    res = res + ...
        sum(((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis)).^2) + ...
        sum(((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis)).^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = compute_res_l1(P,U,u)
res = 0;
for i = 1:length(P);
    %vis = isfinite(u{i}(1,:));
    vis = find(u{i}(1,:)<inf);
    res = res + ...
        sum( abs( ((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis)) ) ) + ...
        sum( abs( ((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis)) ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pnew,Unew] = update_var(d,P,U,u)
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];

dpointvar = [0; d(1:(3*size(U,2)-1))];
dpointvar = reshape(dpointvar, size(U(1:3,:)));
dcamvar = [0;0;0;0;0;0;d(3*size(U,2):end)];
dcamvar = reshape(dcamvar,[6 length(P)]);

Unew = [U(1:3,:) + dpointvar; ones(size(U(1,:)))];

Pnew = cell(size(P));
for i=1:length(P);
    R0 = P{i}(:,1:3);
    t0 = P{i}(:,4);
    R = expm(Ba*dcamvar(1,i) + Bb*dcamvar(2,i) + Bc*dcamvar(3,i))*R0;
    t = t0 + dcamvar(4:6,i);
    Pnew{i} = [R t];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,uStep,UrowIndex] = setup_lin_system(P,U,u)
uStep = zeros(length(u),1);
numpts = size(U,2);
UrowIndex = zeros(numpts,length(u));
%Basis for the tangent plan of the rotation manifold.
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];

da1 = cell(size(u));
db1 = cell(size(u));
dc1 = cell(size(u));
dt11 = cell(size(u));
dt21 = cell(size(u));
dt31 = cell(size(u));
dU11 = cell(size(u));
dU21 = cell(size(u));
dU31 = cell(size(u));
da2 = cell(size(u));
db2 = cell(size(u));
dc2 = cell(size(u));
dt12 = cell(size(u));
dt22 = cell(size(u));
dt32 = cell(size(u));
dU12 = cell(size(u));
dU22 = cell(size(u));
dU32 = cell(size(u));

U = U(1:3,:);

for i=1:length(P);
    %compute the derivatives for both residuals in all images
    %a,b,c - rotation parameters for camera i
    %U1,U2,U3 - 3d point parameters
    %t1,t2,t3 - translation parameters for the camera.
    R0 = P{i}(:,1:3);
    
    t0 = P{i}(:,4);
    da1{i} =  (Ba(1,:)*R0*U)./(R0(3,:)*U+t0(3)) - ...
        (R0(1,:)*U+t0(1))./((R0(3,:)*U+t0(3)).^2).*(Ba(3,:)*R0*U);
    
    da2{i} = (Ba(2,:)*R0*U)./(R0(3,:)*U+t0(3)) - ...
        (R0(2,:)*U+t0(2))./((R0(3,:)*U+t0(3)).^2).*(Ba(3,:)*R0*U);
    
    db1{i} = (Bb(1,:)*R0*U)./(R0(3,:)*U+t0(3)) - ...
        (R0(1,:)*U+t0(1))./((R0(3,:)*U+t0(3)).^2).*(Bb(3,:)*R0*U);
    
    db2{i} = (Bb(2,:)*R0*U)./(R0(3,:)*U+t0(3)) - ...
        (R0(2,:)*U+t0(2))./((R0(3,:)*U+t0(3)).^2).*(Bb(3,:)*R0*U);
    
    dc1{i} = (Bc(1,:)*R0*U)./(R0(3,:)*U+t0(3)) - ...
        (R0(1,:)*U+t0(1))./((R0(3,:)*U+t0(3)).^2).*(Bc(3,:)*R0*U);
    
    dc2{i} = (Bc(2,:)*R0*U)./(R0(3,:)*U+t0(3)) - ...
        (R0(2,:)*U+t0(2))./((R0(3,:)*U+t0(3)).^2).*(Bc(3,:)*R0*U);
    
    dU11{i} = R0(1,1)./(R0(3,:)*U + t0(3)) - ...
        (R0(1,:)*U + t0(1))./((R0(3,:)*U + t0(3)).^2).*R0(3,1);
    
    dU12{i} = R0(2,1)./(R0(3,:)*U + t0(3)) - ...
        (R0(2,:)*U + t0(2))./((R0(3,:)*U + t0(3)).^2).*R0(3,1);
    
    dU21{i} = R0(1,2)./(R0(3,:)*U + t0(3)) - ...
        (R0(1,:)*U + t0(1))./((R0(3,:)*U + t0(3)).^2).*R0(3,2);
    
    dU22{i} = R0(2,2)./(R0(3,:)*U + t0(3)) - ...
        (R0(2,:)*U + t0(2))./((R0(3,:)*U + t0(3)).^2).*R0(3,2);
    
    dU31{i} = R0(1,3)./(R0(3,:)*U + t0(3)) - ...
        (R0(1,:)*U + t0(1))./((R0(3,:)*U + t0(3)).^2).*R0(3,3);
    
    dU32{i} = R0(2,3)./(R0(3,:)*U + t0(3)) - ...
        (R0(2,:)*U + t0(2))./((R0(3,:)*U + t0(3)).^2).*R0(3,3);
    
    dt11{i} = 1./(R0(3,:)*U+t0(3));
    dt12{i} = zeros(size(dt11{i}));
    
    dt21{i} = zeros(size(dt11{i}));
    dt22{i} = 1./(R0(3,:)*U+t0(3));
    
    dt31{i} = -(R0(1,:)*U+t0(1))./((R0(3,:)*U+t0(3)).^2);
    dt32{i} = -(R0(2,:)*U+t0(2))./((R0(3,:)*U+t0(3)).^2);
    
end
row = [];
col = [];
data = [];
resnum = 0;
B = [];
for i = 1:length(u);
    
    vis = find(isfinite(u{i}(1,:)));
    
    uStep(i) = length(vis);
    UrowIndex(vis,i)=resnum+[1:2:2*length(vis)]';
    %The First residual:
    %3D-point parameters:
    %U1-coeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+1]'];
    data = [data; dU11{i}(vis)'];
    %U2-coeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+2]'];
    data = [data; dU21{i}(vis)'];
    %U3-coeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; [vis*3]'];
    data = [data; dU31{i}(vis)'];
    %Camera parameters
    %a-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    data = [data; da1{i}(vis)'];
    %b-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    data = [data; db1{i}(vis)'];
    %c-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    data = [data; dc1{i}(vis)'];
    %t_1-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    data = [data; dt11{i}(vis)'];
    %t_2-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    data = [data; dt21{i}(vis)'];
    %t_3-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    data = [data; dt31{i}(vis)'];
    
    %2nd residual:
    %3D-point parameters:
    %U1-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+1]'];
    data = [data; dU12{i}(vis)'];
    %U2-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+2]'];
    data = [data; dU22{i}(vis)'];
    %U3-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; [vis*3]'];
    data = [data; dU32{i}(vis)'];
    %Camera parameters
    %a-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    data = [data; da2{i}(vis)'];
    %b-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    data = [data; db2{i}(vis)'];
    %c-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    data = [data; dc2{i}(vis)'];
    %t_1-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    data = [data; dt12{i}(vis)'];
    %t_2-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    data = [data; dt22{i}(vis)'];
    %t_3-coeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    data = [data; dt32{i}(vis)'];
    resnum = resnum+2*length(vis);
    
    %Constant termerms
    btmp = zeros(2*length(vis),1);
    %1st residual
    btmp(1:2:end) = (P{i}(1,:)*[U(:,vis); ones(size(U(1,vis)))])./(P{i}(3,:)*[U(:,vis); ones(size(U(1,vis)))])-u{i}(1,vis);
    %2nd residual
    btmp(2:2:end) = (P{i}(2,:)*[U(:,vis); ones(size(U(1,vis)))])./(P{i}(3,:)*[U(:,vis); ones(size(U(1,vis)))])-u{i}(2,vis);
    B = [B; btmp];
end
A = sparse(row,col,data);
%A = zeros(row,col,data);
%Lock the koordinatsystem
%First camera constant and
%first coordinat in the first point constant
A = A(:,[1:3*numpts 3*numpts+7:end]);
A = A(:,[2:end]);
