function delta = compute_update_sparse_l1(A,B,UrowIndex)
[nmpts nmcms]= size(UrowIndex);
%the first cordinate on the first point cordinate , and first carema
%location are fixed
deltaP = zeros((nmcms-1)*6,1);
deltaU = zeros(nmpts*3-1,1);
[A_row_num A_clm_num] = size(A);
res = norm(A*[deltaU;deltaP]-B,1);
gain_last = res;
stop = 0;
iter = 1;
mu  = 1;
maxiter = 20;
J = zeros(A_row_num,(nmcms-1)*6);
while ~stop && iter<maxiter
    iter = iter + 1;
    [deltaU, J] = AllUfromP(A, B,UrowIndex, deltaP);

    b = A*[deltaU;deltaP]-B;
    res = norm(b,1);

    repeat=1;
    while repeat
        fprintf('   Outer LP, iteration # %d \n',repeat)
        % Solve first order approximation in U,
        %d = compute_update(A,-B);
        %d = compute_update_trust(A,-B,mu);
        d = compute_update_interior2(J,-b);
        %d = compute_update_interior_trust(A,-B);
        %d = compute_update_trust2(A,-B, mu);
        %d = compute_update_simplex(A,-B);
        deltaPNew = deltaP+d;
        [deltaUNew J]= AllUfromP(A, B,UrowIndex,deltaPNew);
        b2 = A*[deltaUNew;deltaPNew]-B;
        resnew = norm(b2,1);

        gain= (res-resnew)/gain_last;
        while gain < 0 && norm(d,1) > 1e-6
            d = d/2;
            deltaPNew = deltaPNew+d;
            [deltaUNew J]= AllUfromP(A, B,UrowIndex,deltaPNew);
            mu = mu/2;
            b2 = A*[deltaUNew;deltaPNew]-B;
            resnew = norm(b2,1);
            gain= (res-resnew)/gain_last;
        end
        % Update mu
        n1=0.25;  n2=0.75;  c=2;
        if gain<n1
            mu=n1*(min(norm(d,1),mu));
        end
        if gain>n2 && (abs(norm(d,1)-mu)/mu)<1e-3
            mu=c*mu;
        end

        % Accept update if gain large enough
        fprintf('   gain = %f\n',gain);
        if gain>1e-3
            fprintf('   Gain sufficient, accepting update.\n');
            deltaU = deltaUNew;
            deltaP = deltaPNew;            
            repeat=0;
        else
            fprintf('   Gain not large enough, reducing trust-region. mu value: %f\n', mu);
            repeat=repeat+1;
        end

        % --- Stopping criterion ---
        if mu<1e-6 || abs(res - resnew) < 1e-6;  stop=2; repeat=0;  end % Stop if trust region is really small
    end
    gain_last = res-resnew;
    %while resnew > res
    %if resnew > res
    %    lambda = lambda*2;
        %[Pnew,Unew] = update_var(d,P,U,u);
        %resnew = compute_res(Pnew,Unew,u);
    %end
    %U = Unew;
    %P = Pnew;
    fprintf('%d\t%f\n',iter-1,min(resnew,res));
end
delta=[deltaU; deltaP];

function  [deltaU, J] = AllUfromP(A, B,UrowIndex, deltaP)
[nmpts nmcms]= size(UrowIndex);
%the first cordinate on the first point cordinate , and first carema
%location are fixed
%deltaP = zeros((nmcms-1)*6,1);
deltaU = zeros(nmpts*3-1,1);
[A_row_num A_clm_num] = size(A);

J = zeros(A_row_num,nmcms*6-6);
%first get first U without the first coordinate
indtmp = find(UrowIndex(1,:)>0);
ind = [UrowIndex(1,indtmp) UrowIndex(1,indtmp)+1];
Au = A(ind,1:2);
Ap = A(ind,A_clm_num-nmcms*6+7:A_clm_num);
bu = B(ind,:);
[deltaUtmp basicU Jp] = UfromP(deltaP, Au, Ap, bu);
deltaU(1:2) = deltaUtmp;
J(ind,:)=Jp;
for i = 2:nmpts
    indtmp = find(UrowIndex(i,:)>0);
    ind = [UrowIndex(i,indtmp) UrowIndex(i,indtmp)+1];
    Au = A(ind,3*i-3:3*i-1);
    Ap = A(ind,A_clm_num-nmcms*6+7:A_clm_num);
    bu = B(ind,:);
    [deltaU(3*i-3:3*i-1) basicU Jp] = UfromP(deltaP, Au, Ap, bu);
    J(ind,:)=Jp;
end

function  [deltaU, basic, J_p] = UfromP(deltaP, Au, Ap, bu)

[m,n]=size(Au);
% Form LP
bu = bu - Ap*deltaP;
fv=[sparse(2*n,1); ones(m,1); sparse(2*m,1)];
Av=[ [[-Au Au -speye(m) ]; [Au -Au -speye(m) ]] speye(2*m) ];
bv=[-bu; bu ];

% Construct starting point for simplex phase 2
v=zeros(2*n,1);
res=Av(1:m,1:2*n)*v-bv(1:m);
t=abs(res);
s=[t-res ; t+res];
x0=[v;t;s];
small=1e-6;
basic0=[find((t)>=small )]+2*n;
basic0=[basic0; [find((s)>=2*small )]+m+2*n];
ind=find(t<small);
basic0=[basic0 ; ind+2*n+m ; ind+2*n+m*2];
basic0=sort(basic0);
%[z2,basic2]=simplexfas2b(fv,Av,bv,m,r,n,W,[]);
[z,fval,basic]=simplexphase2(fv,Av,bv,x0,basic0);

% Output
basic=sort(basic); 
deltaU=z(1:n)-z(n+1:2*n); 
num_cam_pars = size(Ap,2);
%compute Jacobian for deltaP
B = Av(:,basic);
Ap_basic = [Ap;-Ap];

Jtmp = inv(B)*Ap_basic;
Jz = zeros(2*n+3*m,num_cam_pars);
Jz(basic,:)=Jtmp;
deltaUdeltaP = Jz(1:n,:)-Jz((n+1):2*n,:);
J_p = Ap + Au*deltaUdeltaP;
1;