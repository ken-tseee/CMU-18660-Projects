clc;clear;
tic;
load('feaSubEImg.mat');
%load('feaSubEOvert.mat');
X=[class{1},class{2}];
[numFeature,numSample]=size(X);
Y=ones(1,numSample);
Y(numSample/2+1:numSample)=(-1)*Y(numSample/2+1:numSample);

[W,C,Lambda,avgAc,stdAc,accuracy] = calculation(X,Y); 
toc;
show_chanWeights(abs(W(:,1)));
