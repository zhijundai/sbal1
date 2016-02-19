function [Pnew] = updateP(P, delta_degree)
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];
Pnew = cell(size(P));
for i=1:length(P);
    R0 = P{i}(:,1:3);
    t0 = P{i}(:,4);
    R = expm(Ba*(rand()-0.5)*delta_degree + Bb*(rand()-0.5)*delta_degree + Bc*(rand()-0.5)*delta_degree)*R0;
    t = t0 + [norm(t0)*(rand()-0.5)*delta_degree;norm(t0)*(rand()-0.5)*delta_degree;norm(t0)*(rand()-0.5)*delta_degree];
    Pnew{i} = [R t];
end