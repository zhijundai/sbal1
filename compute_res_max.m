function res_max = compute_res_max(P,U,u)
res_max = 0;
if size(U,1) == 3
U = [U; ones(1,size(U,2))];
end
for i = 1:length(P);
    vis = isfinite(u{i}(1,:));
    %vis = inlier{i};
    res = max(abs((((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis))) ));
    res2=  max(abs(((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis))));
    res_max = max([res_max res res2]);
end