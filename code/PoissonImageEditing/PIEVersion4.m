% 清理Matlab的环境
close all;
clear;
clc;

% 读入2张图片
SourceImg=imread('images/2/fg.jpg');
SourceImg=im2double(SourceImg);
TargetImg=im2double(imread('images/2/bg.jpg'));
figure,imshow(SourceImg);
uiwait(msgbox({'画个mask图,双击画好的区域结束'}));

% 画mask图
getmask=drawfreehand(gca);%gca检测是否出边界
wait(getmask);
SourceMask=getmask.createMask();
figure,imshow(SourceMask);
% TargetImg = im2double(imread('images/1/bg.jpg'));
% SourceImg = im2double(imread('images/1/fg.jpg'));
% SourceMask = 1-im2bw(imread('images/1/mask.jpg'));

% 设定source将要粘贴在target图中的具体位置，并获取TargetImg的长和宽
position_in_target = [10,150];  %[x,y]
[TargetRows,TargetCols,~] = size(TargetImg);

% b=find(X),X是一个矩阵，查询非零元素的位置，如果X是一个行向量，则返回一个行向量，否则返回一个列向量
[row,col] = find(SourceMask);
[SrcBoundry,~] = bwboundaries(SourceMask,8);

% 计算mask框在source图中的大小
start_pos = [min(col)-1,min(row)-1];
end_pos = [max(col)+1,max(row)+1];
frame_size = end_pos - start_pos;

% 如果position_in_target的位置放置frame将超出Target图的范围，则改变position_in_target,以保证frame不会超出Target图的范围
if(frame_size(1)+position_in_target(1) > TargetCols)
    position_in_target(1) = TargetCols - frame_size(1);
end
if(frame_size(2)+position_in_target(2) > TargetRows)
    position_in_target(2) = TargetRows - frame_size(2);
end

roiTarget=TargetImg(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:);
roiSource=SourceImg(start_pos(2):end_pos(2),start_pos(1):end_pos(1),:);
roiMask = SourceMask(start_pos(2):end_pos(2),start_pos(1):end_pos(1));

[m,n,~] = size(roiMask);

% 获取图像g的梯度场 4个方向
% w1 = [0 -1 0;0 1 0;0 0 0];
% w2 = [0 0 0;0 1 0;0 -1 0];
% w3 = [0 0 0;-1 1 0;0 0 0];
% w4 = [0 0 0;0 1 -1;0 0 0];
w = [-1,1,0];
roiSourceGradx = imfilter(double(roiSource),w);
roiSourceGrady = imfilter(double(roiSource),w');
% roiSourceGradx = roiSourceGradx.*roiMask;
% roiSourceGrady = roiSourceGrady.*roiMask;

% roiSourceGrad1 = imfilter(double(roiSource),w1,'replicate');
% roiSourceGrad2 = imfilter(double(roiSource),w2,'replicate');
% roiSourceGrad3 = imfilter(double(roiSource),w3,'replicate');
% roiSourceGrad4 = imfilter(double(roiSource),w4,'replicate');
% roiSourceGrad1 = roiSourceGrad1.*roiMask;
% roiSourceGrad2 = roiSourceGrad2.*roiMask;
% roiSourceGrad3 = roiSourceGrad3.*roiMask;
% roiSourceGrad4 = roiSourceGrad4.*roiMask;

% 获取背景图片的梯度场
roiTargetGradx = imfilter(double(roiTarget),w);
roiTargetGrady = imfilter(double(roiTarget),w');
% roiTargetGrad1 = imfilter(double(roiSource),w1,'replicate');
% roiTargetGrad2 = imfilter(double(roiSource),w2,'replicate');
% roiTargetGrad3 = imfilter(double(roiSource),w3,'replicate');
% roiTargetGrad4 = imfilter(double(roiSource),w4,'replicate');
% % reverseRoiMask = 1 - roiMask;
% roiTargetGradx = roiTargetGradx.*reverseRoiMask;
% roiTargetGrady = roiTargetGrady.*reverseRoiMask;
% roiTargetGrad1 = roiTargetGrad1.*roiMask;
% roiTargetGrad2 = roiTargetGrad2.*roiMask;
% roiTargetGrad3 = roiTargetGrad3.*roiMask;
% roiTargetGrad4 = roiTargetGrad4.*roiMask;

% 获取融合图像的梯度场
% roiGrad = roiSourceGrad + roiTargetGrad;
% roiGrad1 = roiSourceGrad1 + roiTargetGrad1;
% roiGrad2 = roiSourceGrad2 + roiTargetGrad2;
% roiGrad3 = roiSourceGrad3 + roiTargetGrad3;
% roiGrad4 = roiSourceGrad4 + roiTargetGrad4;
% roiGradx = roiSourceGradx+roiTargetGradx;
% roiGrady = roiSourceGrady+roiTargetGrady;
flag = 1;    % 0:不用混合梯度场    1:用混合梯度场
if flag==0
    for i = 1:m
        for j=1:n
            for k=1:3
                roiGradx(i,j,k) = roiSourceGradx(i,j,k);
                roiGrady(i,j,k) = roiSourceGrady(i,j,k);
              
            end
        end
    end
else
    for i =1:m
        for j=1:n
            for k=1:3
                sourceGrad = roiSourceGradx(i,j,k)^2 + roiSourceGrady(i,j,k)^2;
                targetGrad = roiTargetGradx(i,j,k)^2 + roiTargetGrady(i,j,k)^2;
                if sourceGrad > targetGrad
                    roiGradx(i,j,k) = roiSourceGradx(i,j,k);
                    roiGrady(i,j,k) = roiSourceGrady(i,j,k);
                else
                    roiGradx(i,j,k) = roiTargetGradx(i,j,k);
                    roiGrady(i,j,k) = roiTargetGrady(i,j,k);
                end
            end
        end
    end
end

% for i=1:m
%     for j=1:n
%         for k=1:3
%             if  abs(roiSourceGrad1(i,j,k)) > abs(roiTargetGrad1(i,j,k))
%                 roiGrad1(i,j,k) = roiSourceGrad1(i,j,k);
%             else
%                 roiGrad1(i,j,k) = roiTargetGrad1(i,j,k);
%             end
%             if  abs(roiSourceGrad2(i,j,k)) > abs(roiTargetGrad2(i,j,k))
%                 roiGrad2(i,j,k) = roiSourceGrad2(i,j,k);
%             else
%                 roiGrad2(i,j,k) = roiTargetGrad2(i,j,k);
%             end
%             if  abs(roiSourceGrad3(i,j,k)) > abs(roiTargetGrad3(i,j,k))
%                 roiGrad3(i,j,k) = roiSourceGrad3(i,j,k);
%             else
%                 roiGrad3(i,j,k) = roiTargetGrad3(i,j,k);
%             end
%             if  abs(roiSourceGrad4(i,j,k)) > abs(roiTargetGrad4(i,j,k))
%                 roiGrad4(i,j,k) = roiSourceGrad4(i,j,k);
%             else
%                 roiGrad4(i,j,k) = roiTargetGrad4(i,j,k);
%             end
%         end
%     end
% end


% 获取融合图像的散度
lapx = imfilter(double(roiGradx),w);
lapy = imfilter(double(roiGrady),w');
lap = lapx + lapy;
% lap1 = imfilter(double(roiGrad1),w1,'replicate');
% lap2 = imfilter(double(roiGrad2),w2,'replicate');
% lap3 = imfilter(double(roiGrad3),w3,'replicate');
% lap4 = imfilter(double(roiGrad4),w4,'replicate');
% lap = lap1+lap2+lap3+lap4;

% 构造一个索引矩阵 index_A 
index_A = zeros(m,n);
index = 0;
for i=1:m
    for j=1:n
        if roiMask(i,j) == 1
             index = index + 1;
             index_A(i,j) = index;
        end
    end
end

new_roiTarget = roiTarget;
% 构造A和b
pixels_num = index;
A = spalloc(pixels_num,pixels_num,5*pixels_num);
for i=1:m
    for j=1:n
        if(roiMask(i,j) == 1)
           neighbors = 0;
           if ((i-1)>=1)
                if (roiMask(i-1,j) == 1)
                    neighbors = neighbors + 1;
                end
           end
           if ((j-1)>=1)
                if (roiMask(i,j-1) == 1)
                    neighbors = neighbors + 1;
                end
           end
           if ((i+1)<=m)
                if (roiMask(i+1,j) == 1)
                    neighbors = neighbors + 1;
                end
           end
           if ((j+1)<=n)
                if (roiMask(i,j+1) == 1)
                    neighbors = neighbors + 1;
                end
           end
           % 判断是否是边界点  neighbors==4: roi内部点   neighbors==0: roi外部点
           if (neighbors>0 && neighbors<4)    % 边界点，赋值target对应位置的像素
               ax = index_A(i,j);
               A(ax,ax) = 1; 
               b(ax,1:3) = roiTarget(i,j,:);
           end
           if (neighbors==4)  % 内部点，保持原有散度
               % 构造A
               % div a(i,j) = V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1) - 4V(i,j)
                ax = index_A(i,j);
                ay = index_A(i-1,j);  % V(i-1,j)
                A(ax,ay) = 1;  
                ay = index_A(i+1,j);   % V(i+1,j)
                A(ax,ay) = 1;  
                ay = index_A(i,j-1); % V(i,j-1)
                A(ax,ay) = 1; 
                ay = index_A(i,j+1);  % V(i,j+1)
                A(ax,ay) = 1;
                A(ax,ax) = -4;
                b(ax,1:3) = lap(i,j,:);
                new_roiTarget(i,j,:) = [255 255 255];
           end
        end
    end
end
figure,imshow(new_roiTarget);

RGB = A \ b;   
% RGB = uint8(abs(RGB));   %   结果奇怪???

for i=1:m
    for j=1:n
        ax = index_A(i,j);
        if ax ~= 0
            roiTarget(i,j,:) = RGB(ax,:);
        end
    end
end
figure,imshow(roiTarget);
% roiImage = uint8(abs(roiImage));
resultImage = TargetImg;
% resultImage(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:) = resultImage(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:).*uint8(roiMask);
% roiMask = 1 - roiMask;
% roiImage = roiImage.*double(roiMask);
% resultImage(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:) = resultImage(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:)+roiImage;
resultImage(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:) = roiTarget(:,:,:); %只修改ROI区域内的

figure,imshow(resultImage);

% 保存合成后的图片
imwrite(resultImage,'images/2/result1.jpg');