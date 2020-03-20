clc;
clear all;
close all;

addpath('.\fun\');
load patternMappingriuU4_P24;
mapping2 = getliop('lio');

rootpic = 'G:\Database\Outex_TC_00010\';
SNR0= [0 100 30 15 10 5];

for ii=1:length(SNR0)
    SNR=SNR0(ii);
    picNum = 4320;
    rand('state',0);
    randn('state',0);

    
    for i=1:picNum;
        filename = sprintf('%s\\images\\%06d.ras', rootpic, i-1);
        display(['.... ' num2str(i) ])
        Gray = imread(filename);
        Gray = im2double(Gray);
        if SNR~=0
            Gray=awgn(Gray,10*log10(SNR),'measured'); % Add white Gaussian noise.
        end
        img = (Gray-mean(Gray(:)))/std(Gray(:)); % image normalization, to remove global intensity
        %

            img1=meanmeanfilter(img,3);
            img2=meanmeanfilter(img,5);
% 
            LIOP_3 = LGOP3(img,3,24,mapping2,'x');
            LIOP_5 = LGOP3(img1,5,24,mapping2,'x');
            LIOP_7 = LGOP3(img2,7,24,mapping2,'x');
            

            NLBP_3 = NLBP3(img,3,24,patternMappingriuU4_P24,'x');
            NLBP_5 = NLBP3(img1,5,24,patternMappingriuU4_P24,'x');
            NLBP_7 = NLBP3(img2,7,24,patternMappingriuU4_P24,'x');

            
            LGNBP3 = [LIOP_3 NLBP_3];
            LGNBP5 = [LIOP_5 NLBP_5];
            LGNBP7 = [LIOP_7 NLBP_7];
            

            hist1(i,:) = [LGNBP3 LGNBP5 LGNBP7];

%       

    end

    
    %%
    trainTxt = sprintf('%s000\\train.txt', rootpic);
    testTxt = sprintf('%s000\\test.txt', rootpic);
    [trainIDs, trainClassIDs] = ReadOutexTxt(trainTxt);
    [testIDs, testClassIDs] = ReadOutexTxt(testTxt);

        CP1(ii,:) = cal_AP(hist1,trainIDs, trainClassIDs,testIDs, testClassIDs) %

    
end