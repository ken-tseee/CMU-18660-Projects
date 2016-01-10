function imgOut = imgRecover(imgIn, blkSize, numSample)
% Recover the input image from a small size samples
%
% INPUT:
%   imgIn: input image
%   blkSize: block size
%   numSample: how many samples in each block
%
% OUTPUT:
%   imgOut: recovered image
%
% @ 2011 Huapeng Zhou -- huapengz@andrew.cmu.edu

imgOut = zeros(size(imgIn));
hNum = size(imgIn,1) / blkSize;
wNum = size(imgIn,2) / blkSize;
hB = blkSize;
wB = blkSize;
area = wB * hB;
M = 20;
m = floor(numSample / 6);
lambda = zeros(ceil(numSample / 5),1);
lenLambda = size(lambda,1);
for i=1:size(lambda,1)-1
    lambda(i,1) = lambda(i,1) + i * 5;
end
lambda(lenLambda) = numSample;

for hN=1:hNum
    for wN=1:wNum
        n0 = (hN - 1) * hB + 1;
        nt = hN * hB;
        m0 = (wN - 1) * wB + 1;
        mt = wN * wB;
        blkIn = imgIn(n0:nt,m0:mt);
        T = DCT(blkIn);
        img = reshape(blkIn',area,1);
       
        error = zeros(lenLambda,1);
        for t=1:M
            [B,idxSam] = datasample(img,numSample,'Replace',false);
            A = T(idxSam,:);
            for idxLambda=1:lenLambda
                [testSet,idxTest] = datasample(B,m,'Replace',false);
                
                flag = zeros(size(B,1),1);
                flag(idxTest) = 1;
                flag = ~ flag;
                trainA = A(flag,:);
                trainB = B(flag,:);
                
                trainAlpha = OMP(trainA,trainB,lambda(idxLambda));
                trainC = T * trainAlpha;
                trainSamSet = trainC(idxSam);
                trainSet = trainSamSet(idxTest);
                
                error(idxLambda) = error(idxLambda) + norm(testSet - trainSet);
            end
        end
        error = error ./ M;
        [~,minIdx] = min(error);

        [B,idxSam] = datasample(img,numSample,'Replace',false);
        A = T(idxSam,:);
        alpha = OMP(A,B,lambda(minIdx));
        C = T * alpha;
        lenC = size(C,1);
        blkOut = zeros(hB,wB);
        for i=1:lenC
           xidx = ceil(i / wB);
           yidx = i - (xidx - 1) * wB; 
           blkOut(xidx,yidx) = C(i,1);
        end
        imgOut(n0:nt,m0:mt) = blkOut;
    end
end
lambda(minIdx)
end
