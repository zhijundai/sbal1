function xp = compute_update_interior2(G,y)
%c = G'*y;
%x0 = inv(G'*G)*c;
n=size(G,2);
x0=zeros(n,1);
xp = l1decode_pd2(x0, G, [], y, 1e-6, 10);