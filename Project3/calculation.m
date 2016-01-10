function [W,C,Lambda,avgAc,stdAc,accuracy] = calculation(X,Y)
[numFeature,numSample]=size(X);
class1=X(:,1:numSample/2);
class2=X(:,numSample/2+1:numSample);
label1=Y(1:numSample/2);
label2=Y(numSample/2+1:numSample);

setPara.W=ones(numFeature,1);
setPara.C=0;
setPara.t=1;
setPara.Epsilon=0.000001;
setPara.Beta=15;
setPara.Tmax=1000000;

numFold=6;
accuracy=zeros(numFold,1);
unit=numSample/2/numFold;
Lambda=zeros(numFold,1);
W=zeros(numFeature,numFold);
C=zeros(1,numFold);

for time=1:numFold
    idxStart=(time-1)*unit+1;
    idxEnd=time*unit;
    idxTest=idxStart:idxEnd;
    idxTrain=setdiff(1:numSample/2,idxTest);
    sampleTest=[class1(:,idxTest),class2(:,idxTest)];
    sampleTrain=[class1(:,idxTrain),class2(:,idxTrain)];
	labelTest=[label1(idxTest),label2(idxTest)];
	labelTrain=[label1(idxTrain),label2(idxTrain)];
    
    optLambda=getOptLambda(sampleTrain, labelTrain, setPara);
    
    sizeXi=size(sampleTrain,2);
    Xi=zeros(sizeXi,1);
    for idxXi=1:sizeXi
    	Xi(idxXi,1)=max(0,1-labelTrain(idxXi)*(setPara.W'*sampleTrain(:,idxXi)+setPara.C))+0.001;
    end
    init_Z=[setPara.W;setPara.C;Xi];
    
    t=setPara.t;
    while (t<=setPara.Tmax)
        optZ=solveOptProb_NM(init_Z,optLambda,t,sampleTrain,labelTrain,setPara.Epsilon);
        init_Z=optZ;
        t=setPara.Beta*t;
    end

    optW=optZ(1:numFeature,1);
	optC=optZ(numFeature+1,1);
	labelPredict=optW'*sampleTest+optC;
    predict=labelPredict.*labelTest;
	sizeLabelPridict=length(labelPredict);
    accuracy(time,1)=sum(predict>0)/sizeLabelPridict;
    
    Lambda(time,1)=optLambda;
    W(:,time)=optW;
    C(:,time)=optC;
end

avgAc=mean(accuracy);
stdAc=std(accuracy);

end
