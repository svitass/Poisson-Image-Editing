function [points] = Bresenham(x0,y0,x1,y1)
% 光栅化
dx = x1 - x0; % x偏移量
dy = y1 - y0; % y偏移量
if dx > 0     % x伸展方向
    ux = 1;
else
    ux = -1;
end
if dy > 0     % y伸展方向
    uy = 1;
else
    uy = -1;   
end
dx2 = abs(dx*2);   % x偏移量乘2
dy2 = abs(dy*2);   % y偏移量乘2
points = [[x0,y0]];
if abs(dx) > abs(dy)  % 以x为增量方向计算
    e = -dx;
    x = x0;
    y = y0;
   while x~=x1
        e = e + dy2;
        if e > 0  % e是整数且大于0时表示要取右上的点（否则是右下的点）
            if y~=y1
                y = y + uy;
            end
            e = e - dx2;
        end
        x = x + ux;
        points(end+1,:) = [x,y];
   end
else   % 以y为增量方向计算
    e = -dy;
    x = x0;
    y = y0;
    while y~=y1
        e = e + dx2;
        if e>0
            if x~=x1
                x = x + ux;
            end
            e = e - dy2;
        end
        y = y + uy;
        points(end+1,:)=[x,y];
    end
end
end

