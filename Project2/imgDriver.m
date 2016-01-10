clc;
imgIn = imgRead('fishing_boat.bmp');imgOut = imgRecover(imgIn,8,50);
%imgIn = imgRead('lena.bmp');imgOut = imgRecover(imgIn,16,150);
imgShow(medfilt2(imgOut, [3 3]));
error = mean(mean((imgOut - imgIn) .^ 2));
medianError = mean(mean((medfilt2(imgOut, [3 3]) - imgIn) .^ 2));