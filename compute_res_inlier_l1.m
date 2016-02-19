function res = compute_res_inlier_l1(P,U,u, inlier, display)
if nargin < 5
    display = false;
end
res = 0;
num = 0;
if size(U,1) == 3
U = [U; ones(1,size(U,2))];
end
for i = 1:length(P);
    %vis = isfinite(u{i}(1,:));
    vis = inlier{i};
    num  = num + size(vis,2);
    res = res + ...
        sum( abs( ((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis)) ) ) + ...
        sum( abs( ((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis)) ) );
end
res = res/num;

if display
residue = [];
for i = 1:length(P);
    %vis = isfinite(u{i}(1,:));
    vis = find(u{i}(1,:)<inf);
    %vis = inlier{i};
    %num  = num + size(vis,2);
    residue = [residue abs( ((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis)) ) + ...
            abs( ((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis)) )];
    
end
  figure;
  plot(1:size(residue,2),residue);
end