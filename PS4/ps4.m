%%Q1
load('digitData.mat')
%Q1(a)
m = zeros(size(digitData.train.data,2),1);
for i = 1: size(digitData.train.data,2)
    m(i) = mean(digitData.train.data(:,i));
end
Y = digitData.train.data - ones(size(digitData.train.data,1),1) * m';

%Q1(b)
[Uy, Sy, Vy] = svd(Y,'econ');

%Q1(c)
M = 1/256* (Y'*Y);
[~,D] = eig(M);
total = sum(sum(D));
for i = 256:-1:1
    e(257-i) = D(i,i) / total;
end
figure;plot(e)


%Q1(d)
% figure; imshow(Uy);
% pc = pca(Y);
% data = Uy(:,2) * pc(:,2)';
% data_reshape = reshape(data(1,:), 16, 16);
% figure;imagesc(data_reshape);
% test = Uy * pc(:,1)

%Q1(e)

[~,I] = sort(Vy(:,1));
m1 = zeros(size(digitData.train.data));
for i = 1: size(digitData.train.data,2)
    m1(:,i) = digitData.train.data(:,I(i));
end

m2 = digitData.train.data * Vy(:,1:2);
figure;gscatter(m2(:,1),m2(:,2),digitData.train.digit)
xlabel('v1')
ylabel('v2')


%Q2

