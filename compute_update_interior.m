function xp = compute_update_interior(G,y)
%c = G'*y;
%x0 = inv(G'*G)*c;
n=size(G,2);
x0=zeros(n,1);
xp = l1decode_pd(x0, G, [], y, 1e-6, 30);

% large scale
%gfun = @(z) G*z;
%gtfun = @(z) G'*z;
%xp = l1decode_pd(x0, gfun, gtfun, y, 1e-3, 100, 1e-8, 200);
