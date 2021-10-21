%% Problem Set 
% Chao Gao
% load  data
load('PS7_data.mat');

%% Problem 1
% precompute A^2, A^3
A0 = A^0;
A2 = A^2;
A3 = A^3;
u1 = zeros(size(TF_names));
u2 = zeros(size(TF_names));
trWc = zeros(size(TF_names));
for i = 1:size(TF_names)
    disp(i)
    B_temp = diag(B(:,i));
    C = [B_temp, A*B_temp, A2*B_temp];
    z = X_esc - A3*X_fib(:,8);
    u_star = pinv(C) * z;
    u1(i) = norm(C*u_star - z);
    u2(i) = 3*(transpose(u_star) * u_star);
    trWc(i) = trace(B_temp*transpose(B_temp)*(transpose(A0)*A0+transpose(A)*A+transpose(A2)*A2));
end

control_measure = [u1,u2,trWc];
csvwrite('control_measure.csv',control_measure);

%% Problem 2
[new_order, original_idx] = sort(u1,'ascend');
find(original_idx == 1)
find(original_idx == 2)
find(original_idx == 3)
find(original_idx == 4)
find(original_idx == 5)


[new_order, original_idx] = sort(trWc,'descend');
find(original_idx == 1)
find(original_idx == 2)
find(original_idx == 3)
find(original_idx == 4)
find(original_idx == 5)