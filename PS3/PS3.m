%%Q1
%Q1(a)
A_1 = [1 2; 3 4];
A_2 = [5 6; 7 8];
k = kron(A_1, A_2);
%Q1(b)
t = ttt(tensor(A_1),tensor(A_2));
%Q1(c)
tenmat(t,[3,1],[4,2]);

%%Q2
A = cat(3,A_1,A_2);
T = cp_als(tensor(A),1);
X1 = CP_reconstruct(T.U{1},T.U{2},T.U{3},T.lambda);
T = cp_als(tensor(A),2);
X2 = CP_reconstruct(T.U{1},T.U{2},T.U{3},T.lambda);
T = cp_als(tensor(A),3);
X3 = CP_reconstruct(T.U{1},T.U{2},T.U{3},T.lambda);
T = cp_als(tensor(A),4);
X4 = CP_reconstruct(T.U{1},T.U{2},T.U{3},T.lambda);
error1 = norm(reshape(A-X1,1,[]));
error2 = norm(reshape(A-X2,1,[]));
error3 = norm(reshape(A-X3,1,[]));
error4 = norm(reshape(A-X4,1,[]));

%%Q3

A_mode1 = zeros(2,4);
A_mode2 = zeros(2,4);
A_mode3 = zeros(2,4);
tensor_mode1 = tenmat(A,1);
tensor_mode2 = tenmat(A,2);
tensor_mode3 = tenmat(A,3);
for i = 1:8
    A_mode1(i) = tensor_mode1(i);
    A_mode2(i) = tensor_mode2(i);
    A_mode3(i) = tensor_mode3(i);
end

[~,S1,~] = svd(A_mode1);
[~,S2,~] = svd(A_mode2);
[~,S3,~] = svd(A_mode3);

X_norm = norm(reshape(A,1,[]))^2;
S1_norm = norm(reshape(S1,1,[]))^2;
S2_norm = norm(reshape(S2,1,[]))^2;
S3_norm = norm(reshape(S3,1,[]))^2;

%Proof:
%The norm of tensor X is the same as reshape the tensor into a 2D matrix
%and calculating the matrix's Frobenius norm. The Frobenius norm of the
%matrix is the same with different mode unfolding. The Frobenius norm of a
%2D matrix equals to the sum of the squared singular value because X =
%s1*u1*v1' + s2*u2*v2' and u1,u2,v1,v2 are normalized.

%%Q4
tuck = tucker_als(tensor(A),2);
ho = hosvd(tensor(A),0.01);


x = linspace(-1,1,100);
y = normpdf(x,0,0.5)';
A1 = repmat(y,1,100);
A2 = A1;
A3 = A1;
A2(:,51:100) = A2(:,51:100) + 0.1;
A3(:,51:100) = A3(:,51:100) - 0.1;
X = zeros(90,100,100);
for i = 1:30
X(i,:,:) = A1 + normrnd(0,0.1,[100,100]);
X(i+30,:,:) = A2 + normrnd(0,0.1,[100,100]);
X(i+60,:,:) = A3 + normrnd(0,0.1,[100,100]);
end
% CP decomposition
T = cp_als(tensor(X),1);
figure
plot(T.U{1}) 
% HOSVD
T = tucker_als(tensor(X),[90 100 100]);
figure
plot(T.U{1}(:,1))
