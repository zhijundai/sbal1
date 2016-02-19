function res = compute_res_inlier_l2(P,U,u, inlier, display, max_res, filename)
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
        sum(((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis)).^2) + ...
        sum(((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis)).^2);
end
res = res/num;

if display
residue = [];
isIn = [];
figure;
cur_x = 1;
for i = 1:length(P);
    %vis = isfinite(u{i}(1,:));
    vis = find(u{i}(1,:)<inf);
    in = inlier{i};
    %num  = num + size(vis,2);
    %residue = [residue  sqrt( ( ((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(1,vis)) ).^2 + ...
    %        ( ((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - u{i}(2,vis)) ).^2 )];
    for j = 1:size(vis,2)
        isInlier = true;
        if isempty(find(vis(j)==in))
        isInlier = false;
        end
        res1 = P{i}(1,:)*U(:,vis(j))./(P{i}(3,:)*U(:,vis(j))) - u{i}(1,vis(j));
        res2 =    ( ((P{i}(2,:)*U(:,vis(j)))./(P{i}(3,:)*U(:,vis(j))) - u{i}(2,vis(j))));
        residue = [residue   res1  res2];
        isIn = [isIn isInlier isInlier];
        if isInlier
            line([cur_x cur_x], [0, res1]);
            cur_x = cur_x + 1;
            line([cur_x cur_x], [0, res2]);
        else
            line([cur_x cur_x], [0, res1], 'Color', 'r');
            cur_x = cur_x + 1;
            line([cur_x cur_x], [0, res2], 'Color', 'r');
        end
    end
end

%xlim use max_res
ylim([-max_res max_res]);
%path_img = './';
%set(gcf, 'PaperPosition', [0 0 10 5]);
%print('-depsc',[path_img filename '.eps']);
%print('-dpng',[path_img filename '.png']);
%  figure;
 % plot(1:size(residue,2),residue);
 %x = 1:size(residue,2);
 %line(x(isIn), residue(isIn), 'Color', 'r');
end