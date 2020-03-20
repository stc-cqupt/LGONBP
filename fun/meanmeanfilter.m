function imgCurr=meanmeanfilter(img,H)
filWin = H;
halfWin = (filWin-1)/2;
imgExt = padarray(img,[halfWin halfWin],'symmetric','both');
imgblks = im2col(imgExt,[filWin filWin],'sliding');
%each column of imgblks represent a feature vector
imgMean = mean(imgblks);
imgCurr = reshape(imgMean,size(img));