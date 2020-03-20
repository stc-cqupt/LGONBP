%  ReadOutexTxt gets picture IDs and class IDs from txt for Outex Database
%  [filenames, classIDs] = ReadOutexTxt(txtfile) gets picture IDs and class
%  IDs from TXT file for Outex Database


function [filenames, classIDs] = ReadOutexTxt(txtfile)
% Version 1.0
% Authors: Zhenhua Guo, Lei Zhang and David Zhang
% Copyright @ Biometrics Research Centre, the Hong Kong Polytechnic University

fid = fopen(txtfile,'r');
tline = fgetl(fid); % get the number of image samples
i = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    index = findstr(tline,'.');
    i = i+1;
    filenames(i) = str2num(tline(1:index-1))+1; % the picture ID starts from 0, but the index of Matlab array starts from 1
    classIDs(i) = str2num(tline(index+5:end)); 
end
fclose(fid);