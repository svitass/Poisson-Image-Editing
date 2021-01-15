function [mypath,dis] = dijkstra(a,sb,db,k)
%     a(find(a==0))=inf;  % 将a=0的数全部替换为无穷大(会出现内存不足的错)
    pb(1:length(a))=0;
    d(1:length(a)) = inf;
    prior(1:length(a))=0;
    prior(sb)=sb;
    d(sb)=0;
    tb = find(a(sb,:)~=0);
    tbLen = length(tb);
    for i=1:tbLen    % 更新当前点到其他点的距离
        v = tb(i);
        d(v) = abs(full(a(sb,v))-k);
        prior(v)=sb;
        a(sb,v) = 0;
        a(v,sb) = 0;
    end
    pb(sb)=1; % 当一个点已经求出到原点的最短距离时，其下标i对应的pb(i)赋1
    path = sb;% 存放存入S集合的顺序
    temp = sb;
    length(a)
    for t = 1:length(a)
        % 从未被访问的点中找出最小距离的点
        minn = Inf;
        for i=1:length(a)
            if pb(i)==0 && d(i)<minn
                minn = d(i);
                temp = i;
            end
        end
        pb(temp)=1;
        path = [path,temp];
        tb = find(a(temp,:)~=0);
        tbLen = length(tb);
        for i=1:tbLen
            v = tb(i);
            if d(temp)+abs(full(a(temp,v))-k)< d(v)
                d(v) = d(temp)+abs(full(a(temp,v))-k);
                a(temp,v) = 0;
                a(v,temp) = 0;
                prior(v)=temp;
            end
        end
    end
%     index = find(path==db);
%     mypath = path(1:index);
    v = prior(db);
    mypath=[v,db];
    if v~=0
        while prior(v)~=sb
            v = prior(v);
            mypath=[v,mypath];
        end 
        mypath = [sb,mypath];
    else
        mypath = [];
    end
    dis = d(db);
end
