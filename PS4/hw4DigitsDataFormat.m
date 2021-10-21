% Goal: download and format digits dataset
%
% Resources
%   https://web.stanford.edu/~hastie/ElemStatLearn/datasets/zip.digits/
%   https://umich.instructure.com/courses/274657/assignments/751317
%
% Output
%   digitData:  Digits dataset conatining training and testing data
%               structure with 2 fields (below)
%       train:  Training data. Structure with 2 fields (below)
%           digit:  vector describing the "data" category; (digits, 0-9)
%           data:   matrix of digit data. Data is stored in an N x 256
%                   matrix. Each row corresponds to the category in
%                   "digit". Images can be recapitulated by reshaping the
%                   each row into a 16 x 16 matrix (see example within
%                   script)
%       test:  Testing data. Structure with 2 fields, same as above
%
% Scott Ronquist, scotronq@umich.edu. 4/23/19

%% Set up
clear
close all

% download data
trainOutFile = websave('zip.train.gz',...
    'https://web.stanford.edu/~hastie/ElemStatLearn/datasets/zip.train.gz');
testOutFile = websave('zip.test.gz',...
    'https://web.stanford.edu/~hastie/ElemStatLearn/datasets/zip.test.gz');

% unzip and delete gz
trainOutFileUnzip = gunzip(trainOutFile);
delete(trainOutFile)

testOutFileUnzip = gunzip(testOutFile);
delete(testOutFile)

% read into MATLAB
digitTrain = readtable(trainOutFileUnzip{1},'filetype','text');
delete(trainOutFileUnzip{1})

digitTest = readtable(testOutFileUnzip{1},'filetype','text');
delete(testOutFileUnzip{1})

% format data
digitData.train.digit = digitTrain{:,1};
digitData.train.data = digitTrain{:,2:257};
digitData.test.digit = digitTrain{:,1};
digitData.test.data = digitTrain{:,2:257};

% example - display digits
selectDigit = 3;
digitLoc = find(digitData.train.digit == selectDigit);
figure
for i = 1:16
    subplot(4,4,i)
    dataReshaped = reshape(digitData.train.data(digitLoc(i),:),16,16)';
    imagesc(dataReshaped),
    axis image
end

% save
save('digitData.mat','digitData')
