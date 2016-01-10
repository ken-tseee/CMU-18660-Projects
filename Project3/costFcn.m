function [G, H] = costFcn(Lambda,t,Z,X,Y)
% Compute the cost function F(Z)
%
% INPUTS: 
%   Z: Parameter values
% OUTPUTS
%   F: Function value
%   G: Gradient value
%   H: Hessian value
%
% @ 2011 Kiho Kwak -- kkwak@andrew.cmu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To improve the excution speed, please program your code with matrix
% format. It is 30 times faster than the code using the for-loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X=M*N  N is # of samples, M is # of features
% Y=1*N
% W=M*1
% C=1*1
% Xi=N*1

[numFeature,numSample]=size(X);
lenZ=size(Z,1);
W=Z(1:numFeature,1);
C=Z(numFeature+1,1);
Xi=Z(numFeature+2:lenZ,1);

G=zeros(lenZ,1);
denominator=(W'*X).*Y+C*Y+Xi'-1;

% temp1=sum(Xi)+Lambda*(W'*W);
% temp2=(1/t)*sum(log(denominator));
% temp3=(1/t)*sum(log(Xi));
% F=temp1-temp2-temp3;

G(1:numFeature,1)=2*Lambda*W-(1/t)*X*(Y./denominator)';
G(numFeature+1,1)=-(1/t)*sum((denominator).^(-1).*Y);
G(numFeature+2:lenZ,1)=1-(1/t)*(denominator.^(-1))'-(1/t)*Xi.^(-1);

sY=Y.^2;
sDeno=denominator.^2;

H_WW=zeros(numFeature,numFeature);
for j=1:numFeature
    for k=1:numFeature
       H_WW(j,k)=(1/t)*(X(j,:).*X(k,:))*(sY./sDeno)';
       if (j==k)
          H_WW(j,k)=H_WW(j,k)+2*Lambda; 
       end
    end
end

H_WC=(1/t)*X*(sY./sDeno)';

H_WXi=zeros(numFeature,numSample);
for j=1:numFeature
   for i=1:numSample
       H_WXi(j,i)=(1/t)*X(j,i).*Y(1,i)./sDeno(1,i);
   end
end

H_CC=(1/t)*sum(sY./sDeno);

H_CXi=(1/t)*Y./sDeno; 

XiXi=(1/t)*((1./sDeno)+(Xi'.^(-2)));
H_XiXi=diag(XiXi);

H=[H_WW  , H_WC  , H_WXi ;
   H_WC' , H_CC  , H_CXi ;
   H_WXi', H_CXi', H_XiXi];

end
