load('PS7_data.mat');
%% Problem 1
% precompute A^2, A^3
A_0 = A^0;
A_2 = A^2;
A_3 = A^3;
u1 = zeros(size(TF_names));
u2 = zeros(size(TF_names));
Wc_trace = zeros(size(TF_names));
for i = 1:size(TF_names)
    B_temp = diag(B(:,i));
    C = [B_temp, A*B_temp, A_2*B_temp];
    z = X_esc - A_3*X_fib(:,8);
    u_star = pinv(C) * z;
    u1(i) = norm(C*u_star - z);
    u2(i) = 3*(transpose(u_star) * u_star);
    Wc_trace(i) = trace(B_temp*transpose(B_temp)*(transpose(A_0)*A_0+transpose(A)*A+transpose(A_2)*A_2));
end

control = [u1,u2,Wc_trace];
csvwrite('Problem1.csv',control);
%% Problem 2
[u1_new, u1_idx] = sort(u1,'ascend');
find(u1_idx == 1)
find(u1_idx == 2)
find(u1_idx == 3)
find(u1_idx == 4)
find(u1_idx == 5)


[Wc_new, Wc_idx] = sort(Wc_trace,'descend');
find(Wc_idx == 1)
find(Wc_idx == 2)
find(Wc_idx == 3)
find(Wc_idx == 4)
find(Wc_idx == 5)