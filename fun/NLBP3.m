function result = NLBP3(img,R,P,mapping,mode) %第一个为原始图像，第二个为采样半径，第三个为采样点数，第6个为均值图像

image=img;
d_image=double(image);
radius = R;
neighbors = P;
spoints =zeros(P,2);
% Angle step.
a = 2*pi/P;
for i = 1:neighbors
    spoints(i,1) = -radius*sin((i-1)*a);
    spoints(i,2) = radius*cos((i-1)*a);
end
% Determine the dimensions of the input image.
[ysize xsize] = size(image);



miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));


% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;
% Fill the center pixel matrix C.
C = img(origy:origy+dy,origx:origx+dx);
d_C = double(C);
NLBP_S1=zeros(dy+1,dx+1);
NLBP_S2=zeros(dy+1,dx+1);
NLBP_S3=zeros(dy+1,dx+1);
NLBP_S4=zeros(dy+1,dx+1);
NLBP_C=zeros(dy+1,dx+1);

number = 3;

CC = d_C(:);
[X Y] = size(CC);
CC = sort(CC,'ascend');
x = floor(X*Y/number);

for nur=1:number
    nlm(nur) = mean(CC((nur-1)*(x)+1:nur*(x)));
end

for i = 1:neighbors
    y = spoints(i,1)+origy;
    x = spoints(i,2)+origx;
    % Calculate floors, ceils and rounds for the x and y.
    fy = floor(y); cy = ceil(y); ry = round(y);
    fx = floor(x); cx = ceil(x); rx = round(x);
    % Check if interpolation is needed.
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N = d_image(ry:ry+dy,rx:rx+dx);
        D1{i} = N >= nlm(1);
        D2{i} = N >= nlm(2);
        D3{i} = N >= nlm(3);
%         D4{i} = N >= nlm(4);
        
    else
        % Interpolation needed, use double type images
        ty = y - fy;
        tx = x - fx;
        
        % Calculate the interpolation weights.
        w1 = (1 - tx) * (1 - ty);
        w2 =      tx  * (1 - ty);
        w3 = (1 - tx) *      ty ;
        w4 =      tx  *      ty ;
        % Compute interpolated pixel values
        N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
            w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
        D1{i} = N >= nlm(1);
        D2{i} = N >= nlm(2);
        D3{i} = N >= nlm(3);
%         D4{i} = N >= nlm(4);
        
    end
end

for i=1:neighbors
    % Update the result matrix.
    v = 2^(i-1);
    NLBP_S1 = NLBP_S1 + v*D1{i};
    NLBP_S2 = NLBP_S2 + v*D2{i};
    NLBP_S3 = NLBP_S3 + v*D3{i};
%     NLBP_S4 = NLBP_S4 + v*D4{i};
    
end

NLBP_C = d_C>=mean(d_image(:));

if isstruct(mapping)
    bins = mapping.num;
    sizarray = size(NLBP_S1);
    NLBP_S1 = NLBP_S1(:);
    NLBP_S1 = mapping.table(NLBP_S1+1);
    NLBP_S1 = reshape(NLBP_S1,sizarray);
    
    NLBP_S2 = NLBP_S2(:);
    NLBP_S2 = mapping.table(NLBP_S2+1);
    NLBP_S2 = reshape(NLBP_S2,sizarray);
    
    NLBP_S3 = NLBP_S3(:);
    NLBP_S3 = mapping.table(NLBP_S3+1);
    NLBP_S3 = reshape(NLBP_S3,sizarray);
    
%     NLBP_S4 = NLBP_S4(:);
%     NLBP_S4 = mapping.table(NLBP_S4+1);
%     NLBP_S4 = reshape(NLBP_S4,sizarray);
    
    
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    
    NLBP_S1=hist(NLBP_S1(:),0:(bins-1));
    
    NLBP_S2=hist(NLBP_S2(:),0:(bins-1));
    
    NLBP_S3=hist(NLBP_S3(:),0:(bins-1));
    
%     NLBP_S4=hist(NLBP_S4(:),0:(bins-1));
    
    result = [NLBP_S1 NLBP_S2 NLBP_S3];% NLBP_S4
    
else
    
    idx = find(NLBP_C);
    NLBPC_1 = NLBP_S1;
    NLBPC_1(idx) = NLBPC_1(idx)+bins;
    NLBP_S1=hist(NLBPC_1(:),0:2*bins-1);
    
    NLBPC_2 = NLBP_S2;
    NLBPC_2(idx) = NLBPC_2(idx)+bins;
    NLBP_S2=hist(NLBPC_2(:),0:2*bins-1);
    
    NLBPC_3 = NLBP_S3;
    NLBPC_3(idx) = NLBPC_3(idx)+bins;
    NLBP_S3=hist(NLBPC_3(:),0:2*bins-1);
    
%     NLBPC_4 = NLBP_S4;
%     NLBPC_4(idx) = NLBPC_4(idx)+bins;
%     NLBP_S4=hist(NLBPC_4(:),1:2*bins);

    result = [NLBP_S1 NLBP_S2 NLBP_S3];% NLBP_S4];
end

end