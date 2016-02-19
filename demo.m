%% experiment on Dino data set
if 1 %set to 0 if you do not want to run this
% load data 
load dino.mat

% SBA using the L1 norm
tic
[U_l1,P_l1,ress_l1] = sbal1(U,P,u,50);
toc
% SBA using the L2 norm
tic
[U_l2,P_l2,ress_l2] = bundle(U,P,u,50);
toc
end

%% experiment on the robustness against outliers
if 1 % set to 0 if you do not want to run this
% random data generation
nbrimages=10+ceil(rand()*10);
nbrpoints=100 + ceil(rand()*200);
%noise_level = 0;
noise_level = 0.025;
option={['noise=' int2str(noise_level)],'calibrated','cube=1','distance=2'};
[U,P,u,u0]=randomdata(nbrimages,nbrpoints,option);

% 30 percent of random missing data
for i = 1:length(u)
    num = size(u{i},2);
    for current=1:num        
        if rand() < 0.3
            u{i}(1:3,current) =[inf;inf;inf];
        end
    end
end

running_cases = 1;
res_all_l1 = cell(1,running_cases);
res_all_l2 = cell(1,running_cases);
res_final_inlier =cell(1,running_cases);
U_all_l1 = cell(1,running_cases);
U_all_l2 = cell(1,running_cases);
P_all_l1 = cell(1,running_cases);
P_all_l2 = cell(1,running_cases);
u_outlier_all =  cell(1,running_cases);
times = zeros(running_cases, 2);
for case_id = 1:running_cases
% add noise on U, P
U_in = U;
P_in = P;
u_in = u;
for i= 1:size(U_in,2)
    U_in(1:3,i) = U_in(1:3,i) + 0.05*norm(U_in(1:3,i))*(rand()-0.5);
end
P_in = updateP(P_in,0.05);
% make some measuruments become outliers
inliers=cell(1,length(u));
outliers=cell(1,length(u));

outlier_percentage = 0.05*case_id;
for i = 1:length(u)
    ind = find(u{i}(1,:)<inf);
    num = size(ind,2);
    outlier=[];
    inlier = [];
    for current=1:num        
        if rand() < outlier_percentage
            outlier = [outlier ind(current)];
            % the noise added to u should be much larger than the
            % noise_level
            u_in{i}(1:2,ind(current)) =u{i}(1:2,ind(current)) + (rand(2,1)-0.5)*1;
        else            
            inlier = [inlier ind(current)];
        end
    end
    inliers{i} = inlier;
    outliers{i} = outlier;
    %u{i} = [u{i}(1:2,1:num); ones(1,num)];
end
tic;
[U_l1,P_l1,ress_l1] = sbal1(U_in,P_in,u_in,50);
times(case_id, 1) = toc;
tic;
[U_l2,P_l2,ress_l2] = bundle(U_in,P_in,u_in,50);
times(case_id, 2) = toc;
max_res = compute_res_max(P_l1,U_l1,u_in);
max_res2 = compute_res_max(P_l2,U_l2,u_in);
max_res = max(max_res,max_res2);
displayFigure = true;
res11 = compute_res_inlier_l1(P_l1,U_l1,u_in, inliers, false);
res12 = compute_res_inlier_l2(P_l1,U_l1,u_in, inliers, displayFigure,max_res, sprintf('fig-outlier-l1-%02d', outlier_percentage*100));
res21 = compute_res_inlier_l1(P_l2,U_l2,u_in, inliers, false);
res22 = compute_res_inlier_l2(P_l2,U_l2,u_in, inliers, displayFigure,max_res, sprintf('fig-outlier-l2-%02d', outlier_percentage*100));
res_all_l1{1,case_id} = ress_l1;
res_all_l2{1,case_id} = ress_l2;
res_final_inlier{1,case_id} = [res11, res12; res21, res22];
end
end