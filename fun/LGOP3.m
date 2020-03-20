function  lgop= LGOP3(img,R,P,mapping,mode)

%分成6组
% Check number of input arguments.
d_image=double(img);

radius = R;
neighbors = P;

spoints =zeros(P,2);

% Angle step.
a1 = 2*pi/P;

for i= 1:P
    spoints(i,1) = -radius*sin((i-1)*a1);
    spoints(i,2) = radius*cos((i-1)*a1);
end

neighbors = size(spoints,1);

% Determine the dimensions of the input image.
[ysize xsize] = size(img);


%找出采样点横坐标和纵坐标中的最大值和最小值
miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;  %R=4,bsizey=9
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;  %R=4,bsizey=9

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));   %R=4,origy=5
origx=1-floor(min(minx,0));   %R=4,origx=5

% Minimum allowed size for the input image depends

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = img(origy:origy+dy,origx:origx+dx);
d_C = double(C);


% Initialize the result matrix with zeros.

LGOP_C=zeros(dy+1,dx+1);


%先求领域像素的梯度

for i = 1:P
    y = spoints(i,1)+origy;
    x = spoints(i,2)+origx;
    % Calculate floors, ceils and rounds for the x and y.
    fy = floor(y); cy = ceil(y); ry = round(y);
    fx = floor(x); cx = ceil(x); rx = round(x);
    % Check if interpolation is needed.
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N = d_image(ry:ry+dy,rx:rx+dx);
        N1(i,:)=N(:)';
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
        N1(i,:)=N(:)';
    end
end

N1 = N1';
DN1 = abs(N1 - repmat(mean(N1,2),1,size(N1,2)));
[max_a, index] = max(DN1');

for i = 1:length(max_a)
    
    N1(i,:)= circshift(N1(i,:),-index(i)+1,2);
%     a = N1(i,:);
%     N1(i,:) = [a(index(i):end) a(1:index(i)-1)];
end
N1 = N1';

M1 = N1;
M6 = [M1(6,:);M1(12,:);M1(18,:);M1(24,:)];
M5 = [M1(5,:);M1(11,:);M1(17,:);M1(23,:)];
M4 = [M1(4,:);M1(10,:);M1(16,:);M1(22,:)];
M3 = [M1(3,:);M1(9,:);M1(15,:);M1(21,:)];
M2 = [M1(2,:);M1(8,:);M1(14,:);M1(20,:)];
M1 = [M1(1,:);M1(7,:);M1(13,:);M1(19,:)];
[~, n] = size(M1);
M = [M1 M2 M3 M4 M5 M6];
[~,M] = sort(M,'descend' );
D1_nur = M(:,1:n)';
D2_nur = M(:,n*1+1:n*2)';
D3_nur = M(:,n*2+1:n*3)';
D4_nur = M(:,n*3+1:n*4)';
D5_nur = M(:,n*4+1:n*5)';
D6_nur = M(:,n*5+1:n*6)';

% [~,D1]=sort(M1,'descend' );
% D1_nur = D1';
% [~,D2]=sort(M2,'descend' );
% D2_nur = D2';
% [~,D3]=sort(M3,'descend' );
% D3_nur = D3';
% [~,D4]=sort(M4,'descend' );
% D4_nur = D4';
% [~,D5]=sort(M5,'descend' );
% D5_nur = D5';
% [~,D6]=sort(M6,'descend' );
% D6_nur = D6';

weight = [1000 100 10 1]';
LIOP_1 = D1_nur*weight;
LIOP_2 = D2_nur*weight;
LIOP_3 = D3_nur*weight;
LIOP_4 = D4_nur*weight;
LIOP_5 = D5_nur*weight;
LIOP_6 = D6_nur*weight;


% BRINT_C
LGOP_C = d_C>=mean(d_image(:));


%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    
    sizarray = size(LGOP_C);
    LIOP_1 = mapping.table(LIOP_1);
    LIOP_1 = reshape(LIOP_1,sizarray);
    LIOP_2 = mapping.table(LIOP_2);
    LIOP_2 = reshape(LIOP_2,sizarray);
    LIOP_3 = mapping.table(LIOP_3);
    LIOP_3 = reshape(LIOP_3,sizarray);
    LIOP_4 = mapping.table(LIOP_4);
    LIOP_4 = reshape(LIOP_4,sizarray);
    LIOP_5 = mapping.table(LIOP_5);
    LIOP_5 = reshape(LIOP_5,sizarray);
    LIOP_6 = mapping.table(LIOP_6);
    LIOP_6 = reshape(LIOP_6,sizarray);
end


if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    
    LIOP_1= hist(LIOP_1(:),0:bins-1);
    LIOP_2= hist(LIOP_2(:),0:bins-1);
    LIOP_3= hist(LIOP_3(:),0:bins-1);
    LIOP_4= hist(LIOP_4(:),0:bins-1);
    LIOP_5= hist(LIOP_5(:),0:bins-1);
    LIOP_6= hist(LIOP_6(:),0:bins-1);
    
    
else
    idx = find(LGOP_C);
    LIOPC_1 = LIOP_1;
    LIOPC_1(idx) = LIOPC_1(idx)+bins;
    LIOP_1=hist(LIOPC_1(:),0:2*bins-1);
    
    LIOPC_2 = LIOP_2;
    LIOPC_2(idx) = LIOPC_2(idx)+bins;
    LIOP_2=hist(LIOPC_2(:),0:2*bins-1);
    
    LIOPC_3 = LIOP_3;
    LIOPC_3(idx) = LIOPC_3(idx)+bins;
    LIOP_3=hist(LIOPC_3(:),0:2*bins-1);
    
    LIOPC_4 = LIOP_4;
    LIOPC_4(idx) = LIOPC_4(idx)+bins;
    LIOP_4=hist(LIOPC_4(:),0:2*bins-1);
    
    LIOPC_5 = LIOP_5;
    LIOPC_5(idx) = LIOPC_5(idx)+bins;
    LIOP_5=hist(LIOPC_5(:),0:2*bins-1);
    
    LIOPC_6 = LIOP_6;
    LIOPC_6(idx) = LIOPC_6(idx)+bins;
    LIOP_6=hist(LIOPC_6(:),0:2*bins-1);
    
end

lgop = [LIOP_1 LIOP_2 LIOP_3 LIOP_4 LIOP_5 LIOP_6];

% [max_lgop1, index1] = max(LIOP_1);
% % LIOP_1 = [LIOP_1(index1:end) LIOP_1(1:index1-1)];
% [max_lgop2, index2] = max(LIOP_2);
% % LIOP_2 = [LIOP_2(index2:end) LIOP_2(1:index2-1)];
% [max_lgop3, index3] = max(LIOP_3);
% % LIOP_3 = [LIOP_3(index3:end) LIOP_3(1:index3-1)];
% [max_lgop4, index4] = max(LIOP_4);
% % LIOP_4 = [LIOP_4(index4:end) LIOP_4(1:index4-1)];
% [max_lgop5, index5] = max(LIOP_5);
% % LIOP_5 = [LIOP_5(index5:end) LIOP_5(1:index5-1)];
% [max_lgop6, index6] = max(LIOP_6);
% % LIOP_6 = [LIOP_6(index6:end) LIOP_6(1:index6-1)];
% 
% [~, index] = max([max_lgop1 max_lgop2 max_lgop3 max_lgop4 max_lgop5 max_lgop6]);
% 
% if index == 1
%     lgop = [LIOP_1 LIOP_2 LIOP_3 LIOP_4 LIOP_5 LIOP_6];
% elseif index == 2
%     lgop = [LIOP_2 LIOP_3 LIOP_4 LIOP_5 LIOP_6 LIOP_1];
% elseif index == 3
%     lgop = [LIOP_3 LIOP_4 LIOP_5 LIOP_6 LIOP_1 LIOP_2];
% elseif index == 4
%     lgop = [LIOP_4 LIOP_5 LIOP_6 LIOP_1 LIOP_2 LIOP_3];
% elseif index == 5
%     lgop = [LIOP_5 LIOP_6 LIOP_1 LIOP_2 LIOP_3 LIOP_4];
% else 
%     lgop = [LIOP_6 LIOP_1 LIOP_2 LIOP_3 LIOP_4 LIOP_5];
% end
