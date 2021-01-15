% %% GrabCut 函数应用
% % 从图像中抠出前景色
% 
% clear;clc;
% %读入RGB图像
% SourceImg = imread('images/1/fg.jpg');%可以自行更改图片
% TargetImg=im2double(imread('images/1/bg.jpg'));
% %超像素分割——创建一个同样大小的图像
% L = superpixels(SourceImg,500);
% 
% %% 创建可变多边形
% imshow(SourceImg)
% %①手动创建，gca为函数句柄
% % 在弹出figure窗口直接鼠标点击选取闭环区域。
% h1 = impoly(gca); %impoly() = drawpolygon() = drawpolyline()
% 
% %②定点创建，gca为函数句柄
% % h1 = impoly(gca,[72,105; 1,231; 0,366; 104,359;...
% %         394,307; 518,343; 510,39; 149,72]);
% 
% %% 获取ROI位置坐标
% ROIPoints = getPosition(h1);
% 
% %% 将感兴趣区域（ROI）多边形转换为区域蒙版
% % poly2mask（x，y，m，n）从在x和y处具有顶点的ROI多边形计算大小为m×n的二进制感兴趣区域（ROI）掩模BW。 
% ROI = poly2mask(ROIPoints(:,1),ROIPoints(:,2),size(L,1),size(L,2));
% 
% %% 开始分割
% object = grabcut(SourceImg,L,ROI);
% figure,imshow(object);
% title('物体的边界');

%% Drag and Drop Pasting  内边界  
clear;clc;
%读入RGB图像
SourceImg = imread('images/1/fg.jpg');%可以自行更改图片
TargetImg=im2double(imread('images/1/bg.jpg'));
figure,imshow(SourceImg);
uiwait(msgbox({'内边界：画个mask图,双击画好的区域结束'}));

% 画mask图
getmask=drawfreehand(gca);%gca检测是否出边界
wait(getmask);
object=getmask.createMask();
figure,imshow(object);
title('object的边界');

%% Drag and Drop Pasting  
figure,imshow(SourceImg);
uiwait(msgbox({'画个mask图,双击画好的区域结束'}));

% 画mask图
getmask=drawfreehand(gca);%gca检测是否出边界
wait(getmask);
SourceMask=getmask.createMask();
figure,imshow(SourceMask);
title('ROI区域的边界');

% 获取物体边界和ROI边界围成的区域
band = xor(object,SourceMask);
figure,imshow(band);
title('band');

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

SourceImg = im2double(SourceImg);
TargetImg = im2double(TargetImg);
SourceMask = im2double(SourceMask);
roiTarget=TargetImg(position_in_target(2):position_in_target(2)+frame_size(2),position_in_target(1):position_in_target(1)+frame_size(1),:);
roiSource=SourceImg(start_pos(2):end_pos(2),start_pos(1):end_pos(1),:);
roiMask = SourceMask(start_pos(2):end_pos(2),start_pos(1):end_pos(1));
roiBandMask = band(start_pos(2):end_pos(2),start_pos(1):end_pos(1));
objectMask = object(start_pos(2):end_pos(2),start_pos(1):end_pos(1));


[m,n,~] = size(roiBandMask);

% 构造索引矩阵index_A ，存储band中的像素在稀疏矩阵中的索引
index_A = zeros(m,n);
index = 0;
for i=1:m
    for j=1:n
        if roiBandMask(i,j) == 1
             index = index + 1;
             index_A(i,j) = index;
        end
    end
end

% 构造图，图中的顶点：band中的像素，边的权重：指向像素在源图片和目标图片的像素差，两个不能到达的顶点的边的表示：0
pixels_num = index;
graph = spalloc(pixels_num,pixels_num,5*pixels_num);
for i=1:m
    for j=1:n
        if(roiBandMask(i,j) == 1)
            ax = index_A(i,j);
            ft = roiTarget(i,j,:);
            fs = roiSource(i,j,:);
            w = ft - fs;
            %  (f(p)− f(q)) iscomputedastheL2-norm in color spaces.
%            w = sqrt((w(1,1,1))^2+(w(1,1,2))^2+(w(1,1,3))^2);
%            graph(ax,ax) = w;
            graph(ax,ax) = 0;
           if ((i-1)>=1)
                if (roiBandMask(i-1,j) == 1)
                    ft = roiTarget(i-1,j,:);
                    fs = roiSource(i-1,j,:);
                    w = ft - fs;
                    w = sqrt((w(1,1,1))^2+(w(1,1,2))^2+(w(1,1,3))^2)+1;
                    ay = index_A(i-1,j);
                    graph(ax,ay) = w;
                end
           end
           if ((j-1)>=1)
                if (roiBandMask(i,j-1) == 1)
                    ft = roiTarget(i,j-1,:);
                    fs = roiSource(i,j-1,:);
                    w = ft - fs;
                    w = sqrt((w(1,1,1))^2+(w(1,1,2))^2+(w(1,1,3))^2)+1;
                    ay = index_A(i,j-1);
                    graph(ax,ay) = w;
                end
           end
           if ((i+1)<=m)
                if (roiBandMask(i+1,j) == 1)
                    ft = roiTarget(i+1,j,:);
                    fs = roiSource(i+1,j,:);
                    w = ft - fs;
                    w = sqrt((w(1,1,1))^2+(w(1,1,2))^2+(w(1,1,3))^2)+1;
                    ay = index_A(i+1,j);
                    graph(ax,ay) = w;
                end
           end
           if ((j+1)<=n)
                if (roiBandMask(i,j+1) == 1)
                    ft = roiTarget(i,j+1,:);
                    fs = roiSource(i,j+1,:);
                    w = ft - fs;
                    w = sqrt((w(1,1,1))^2+(w(1,1,2))^2+(w(1,1,3))^2)+1;
                    ay = index_A(i,j+1);
                    graph(ax,ay) = w;
                end
           end
        end
    end
end

% 保存边界点的像素
RoiBoundary = bwboundaries(roiMask,8);   % roi的边界
ObjBoundary = bwboundaries(objectMask,8);   % 物体的边界
RoiBoundary = RoiBoundary{1};
ObjBoundary = ObjBoundary{1};

roiLen = length(RoiBoundary);
objLen = length(ObjBoundary);

%%测试边界
for i=1:roiLen
    pos = RoiBoundary(i,:);
    roiSource(pos(1,1),pos(1,2),:) = [0,0,0];
end
for i=1:objLen
    pos = ObjBoundary(i,:);
    roiSource(pos(1,1),pos(1,2),:) = [0,0,0];
end


% 计算线段C:连接内边界到外边界的最短直线（欧式距离）
pixel_roi = RoiBoundary(1,:);
pixel_obj = ObjBoundary(1,:);
minDist = sqrt(sum((pixel_roi - pixel_obj).^2));
for i = 1:roiLen
    roiPixel = RoiBoundary(i,:);
    for j = 1:objLen
        objPixel = ObjBoundary(j,:);
        dist = sqrt(sum((roiPixel - objPixel).^2));
        if dist < minDist
            minDist = dist;
            pixel_roi = roiPixel;
            pixel_obj = objPixel;
        end
    end
end

% 将线段C光栅化为2D图像的像素列表
points = Bresenham(pixel_roi(1,1),pixel_roi(1,2),pixel_obj(1,1),pixel_obj(1,2));
pointsLen = length(points);

%%展示光栅化结果
% for i=1:pointsLen
%     pos = points(i,:);
%     roiSource(pos(1,1),pos(1,2),:) = [255,0,0];
% end
% figure,imshow(roiSource);
% title("展示内外边界+光栅化");

% 获取邻接点：points作为黄色的点，points下方的点（i+1）或者左侧的点(j-1)作为蓝色的点
pixels_side = zeros(m,n);    % 用于标记蓝色的点
y = abs(pixel_roi(1,1) - pixel_obj(1,1));
x = abs(pixel_roi(1,2) - pixel_obj(1,2));
if y/x <= 1   % 小于pi/4，i+1
    disp('i+1');
    points_blue = Bresenham(pixel_roi(1,1)+1,pixel_roi(1,2),pixel_obj(1,1)+1,pixel_obj(1,2));
    for i=1:length(points_blue)
        ax = points_blue(i,1);
        ay = points_blue(i,2);
        pixels_side(ax,ay) = 1;
        roiSource(ax,ay,:) = [0,0,255];
    end
    for i = 1:pointsLen
        point = points(i,:);
        roiSource(point(1,1),point(1,2),:) = [255,255,0];
    end

else   % j-1
    disp('j-1');
    points_blue = Bresenham(pixel_roi(1,1),pixel_roi(1,2)-1,pixel_obj(1,1),pixel_obj(1,2)-1);
    for i=1:length(points_blue)
        ax = points_blue(i,1);
        ay = points_blue(i,2);
        pixels_side(ax,ay) = 1;
        roiSource(ax,ay,:) = [0,0,255];
    end
    for i = 1:pointsLen
        point = points(i,:);
        roiSource(point(1,1),point(1,2),:) = [255,255,0];
    end
end


% 切断线段两侧点的连通性
adjPoints = [];    % 每一行：[黄点，蓝点]  （graph中的索引）
for i = 1:pointsLen
    point = points(i,:);
    ax = index_A(point(1,1),point(1,2));
    if ax~=0   % 可能光栅化的点不在band内
        if (pixels_side(point(1,1)-1,point(1,2)) == 1)   % 蓝色的点，切断连通性
            ay = index_A(point(1,1)-1,point(1,2));
            if ay~=0
                graph(ax,ay) = 0;
                adjPoints(end+1,:) = [ax,ay]; 
            end
        end
        if (pixels_side(point(1,1)+1,point(1,2)) == 1)
            ay = index_A(point(1,1)+1,point(1,2));
            if ay~=0
                graph(ax,ay) = 0;
                adjPoints(end+1,:) = [ax,ay]; 
            end
        end
        if (pixels_side(point(1,1),point(1,2)-1) == 1)
             ay = index_A(point(1,1),point(1,2)-1);
             if ay~=0
                graph(ax,ay) = 0;
                adjPoints(end+1,:) = [ax,ay]; 
             end
        end
        if (pixels_side(point(1,1),point(1,2)+1) == 1)
             ay = index_A(point(1,1),point(1,2)+1);
             if ay~=0
                graph(ax,ay) = 0;
                adjPoints(end+1,:) = [ax,ay];
             end
        end
    end
end
figure,imshow(roiSource);
title("展示内外边界+光栅化+黄色蓝色点");

% 初始化contour为roi的边界，并且计算k
wSum = 0;
path = [];
for i = 1:roiLen
    boundary = RoiBoundary(i,:);
    ax = boundary(1,1);
    ay = boundary(1,2);
    ft = roiTarget(ax,ay,:);
    fs = roiSource(ax,ay,:);
    w = ft - fs;
    w = sqrt((w(1,1,1))^2+(w(1,1,2))^2+(w(1,1,3))^2);
    wSum = wSum + w;
    vertex = index_A(ax,ay);
    path(end+1) = vertex;
end
k = wSum / roiLen;

% 计算roi边界的路径和 需要用到k
pathLen = length(path);
s = path(pathLen);
e = path(1);
w = full(graph(s,e));  % full:稀疏矩阵转换为全元素矩阵
lossSum = abs(w-k);
for i=1:pathLen-1
    s = path(i);
    e = path(i+1);
    lossSum = lossSum + abs(w-k);
end

% 计算adjPoints中点与点之间的最短路径 
tic
adjPointsLen = length(adjPoints);
minLoss = lossSum;
minPath = path;
fprintf('要循环 %d 次Dijistra\n',adjPointsLen);
i = 1;
while adjPointsLen~=0
    startVertex = adjPoints(1,1);
    endVertex = adjPoints(1,2);
    fprintf('第 %d 次Dijistra\n',i);
    [path,dis] = dijkstra(graph,startVertex,endVertex,k);
    w = full(graph(endVertex,startVertex));
    dis = dis+abs(w-k);
    minLoss
    dis
    if minLoss > dis && (length(path)~=0)
        disp('更新路径');
        minLoss = dis;
        minPath = path;
    end
    % 更新k
    pathLen = length(path);
    if pathLen ~= 0
       s = path(pathLen);
       e = path(1);
       wSum = graph(s,e);
       for j=1:pathLen-1
           s = path(j);
           e = path(j+1);
           w = graph(s,e);
           wSum = wSum + w;
       end
       k = wSum / pathLen; 
    end
    adjPoints = adjPoints(2:end,:);
    adjPointsLen = length(adjPoints);
    i = i+1;
end
toc
%显示优化后的边界
pathLen = length(minPath);
 for j=1:pathLen
    s = minPath(j);
    [x,y] = find(index_A == s);
    roiSource(x,y,:) = [255,0,0];
 end
figure,imshow(roiSource);
title("展示内外边界+光栅化+黄色蓝色点,优化后的边界");

% 根据边界minPath生成Mask：深度优先遍历
roiMask = zeros(m,n);
pathLen = length(minPath);
for j=1:pathLen
    s = minPath(j);
    [x,y] = find(index_A == s);
    roiMask(x,y)=1;
end
queue = pixel_obj;
vis = zeros(m,n);
while length(queue) ~= 0 
   point = queue(end-1:end);
   vis(point(1,1),point(1,2)) = 1;
   queue = queue(1:end-2);
   flag = 0;   % point不是边界点，可向queue中添加邻接点
      vertex = index_A(point(1,1),point(1,2));
      if vertex~=0
        for j=1:pathLen
            s = minPath(j);
            if (vertex == s)
                flag = 1;
                break;
            end
         end
      end
   
      if flag==0
            roiMask(point(1,1),point(1,2)) = 1;
            point1 = [point(1,1),point(1,2)-1];
            if vis(point1(1,1),point1(1,2)) == 0
                queue = [queue,point1];
            end
            point2 = [point(1,1),point(1,2)+1];
            if vis(point2(1,1),point2(1,2)) == 0
                queue = [queue,point2];
            end
            point3 = [point(1,1)-1,point(1,2)];
            if vis(point3(1,1),point3(1,2)) == 0
                queue = [queue,point3];
            end
            point4 = [point(1,1)+1,point(1,2)];
            if vis(point4(1,1),point4(1,2)) == 0
                queue = [queue,point4];
            end
      end 
end
figure,imshow(roiMask);
title("优化后的边界构成的Mask");


%% PossionImageEditing
roiSource=SourceImg(start_pos(2):end_pos(2),start_pos(1):end_pos(1),:);
% 获取图像g的梯度场 4个方向
w = [-1,1,0];
roiSourceGradx = imfilter(double(roiSource),w);
roiSourceGrady = imfilter(double(roiSource),w');
% 获取背景图片的梯度场
roiTargetGradx = imfilter(double(roiTarget),w);
roiTargetGrady = imfilter(double(roiTarget),w');

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

% 获取融合图像的散度
lapx = imfilter(double(roiGradx),w);
lapy = imfilter(double(roiGrady),w');
lap = lapx + lapy;

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
imwrite(resultImage,'images/1/result2.jpg');