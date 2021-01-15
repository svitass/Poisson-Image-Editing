function [mypath,dis] = dijkstra(a,sb,db,k)
%     a(find(a==0))=inf;  % ��a=0����ȫ���滻Ϊ�����(������ڴ治��Ĵ�)
    pb(1:length(a))=0;
    d(1:length(a)) = inf;
    prior(1:length(a))=0;
    prior(sb)=sb;
    d(sb)=0;
    tb = find(a(sb,:)~=0);
    tbLen = length(tb);
    for i=1:tbLen    % ���µ�ǰ�㵽������ľ���
        v = tb(i);
        d(v) = abs(full(a(sb,v))-k);
        prior(v)=sb;
        a(sb,v) = 0;
        a(v,sb) = 0;
    end
    pb(sb)=1; % ��һ�����Ѿ������ԭ�����̾���ʱ�����±�i��Ӧ��pb(i)��1
    path = sb;% ��Ŵ���S���ϵ�˳��
    temp = sb;
    length(a)
    for t = 1:length(a)
        % ��δ�����ʵĵ����ҳ���С����ĵ�
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
