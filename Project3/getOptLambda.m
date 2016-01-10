function [optLambda] = getOptLambda(X, Y, setPara)
% Get the optimal lamda
%
% INPUTS:
%   X(MxN) : trData(i,j) is the i-th feature from the j-th trial
%   Y(Nx1): trData(j) is the label of the j-th trial (1 or -1)
%   setPara : Initialized parameters
%            setPara.t      
%            setPara.zeta   
%            setPara.Tmax   
%            setPara.tol    
%            setPara.W      
%            setPara.C      
%
% OUTPUTS:
%   optiLamda: Optimal lamda value 
%
% @ 2011 Kiho Kwak -- kkwak@andrew.cmu.edu

[numFeature,numSample]=size(X);
W=setPara.W;
C=setPara.C;
Beta=setPara.Beta;
Epsilon=setPara.Epsilon;
Tmax=setPara.Tmax;
optLambda=0;

optional_Lambda = [0.01,1,100,10000];
lenLambda=length(optional_Lambda);
numFold=5;

class1=X(:,1:numSample/2);
class2=X(:,numSample/2+1:numSample);
label1=Y(1:numSample/2);
label2=Y(numSample/2+1:numSample);
avgAccuracy=zeros(lenLambda,1);
unit=numSample/2/numFold;

for idxLambda=1:lenLambda
    Lambda=optional_Lambda(idxLambda);
    accuracy=zeros(numFold,1);
    maxAccuracy=0;
    for time=1:numFold
        idxStart=(time-1)*unit+1;
        idxEnd=time*unit;
        idxTest=idxStart:idxEnd;
        idxTrain=setdiff(1:numSample/2,idxTest);
        sampleTest=[class1(:,idxTest),class2(:,idxTest)];
        sampleTrain=[class1(:,idxTrain),class2(:,idxTrain)];
        labelTest=[label1(idxTest),label2(idxTest)];
        labelTrain=[label1(idxTrain),label2(idxTrain)];
        
        sizeXi=size(sampleTrain,2);
        Xi=zeros(sizeXi,1);
        for idxXi=1:sizeXi
            Xi(idxXi,1)=max(0,1-labelTrain(idxXi)*(W'*sampleTrain(:,idxXi)+C))+0.001;
        end
        init_Z=[W;C;Xi];
        t=setPara.t;
        while (t<=Tmax)
            Z_new=solveOptProb_NM(init_Z,Lambda,t,sampleTrain,labelTrain,Epsilon);
            init_Z=Z_new;
            t=Beta*t;
        end
        W_new=Z_new(1:numFeature,1);
        C_new=Z_new(numFeature+1,1);
        labelPredict=W_new'*sampleTest+C_new;
        predict=labelPredict.*labelTest;
        sizeLabelPridict=length(labelPredict);
        accuracy(time,1)=sum(predict>0)/sizeLabelPridict;
    end
    avgAccuracy(idxLambda,1)=mean(accuracy);
    if (avgAccuracy(idxLambda,1)>maxAccuracy)
        maxAccuracy=avgAccuracy(idxLambda,1);
        optLambda=Lambda;
    end
end

end
