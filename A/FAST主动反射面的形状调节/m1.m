
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=300;            %基准态反射面半径
D=500;            %口径为D时的球面（基准球面）
d=300;            %工作抛物面口径
F=0.466*R;        %焦面 与基准球面半径差
dk=1;             %馈源舱接收信号有效区域直径
alpha=36.795/180*pi;          %观测天体的alpha角
beta=78.169/180*pi;           %观测天体的beta角
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=R-F;            %焦面半径
P=[-r*cos(alpha)*cos(beta),-r*sin(alpha)*cos(beta),-r*sin(beta)];%馈源舱圆形坐标
O_d=P*(R^2-d^2/4)^0.5/r;%工作抛物面口径圆心坐标

l=[];
for i=1:size(points0,1)
    l=[l;norm(points0(i,2:4)-O_d)];
end
points0=sortrows([points0 l],17);
points0=points0(:,1:16);

points_work=[];%工作抛物面内的点
for i=1:size(points0,1)
    if norm(O_d-points0(i,2:4))<d/2
        points_work=[points_work;points0(i,1)];
    end
end

ls_work=lines0;%工作抛物面内的线
for i=size(ls_work,1):-1:1
    if ~ ismember(ls_work(i,2),points_work) || ~ismember(ls_work(i,6),points_work)
        ls_work(i,:)=[];
    end
end
ls_work=ls_work(:,1);

ls_notwork=lines0;%非工作抛物面内的线
for i=size(ls_work,1):-1:1
    ls_notwork(ls_work(i),:)=[];
end
ls_notwork=ls_notwork(:,1);

sfs_work=surfaces0;%工作抛物面内的面
for i=size(sfs_work,1):-1:1
    if ~ismember(sfs_work(i,2),points_work) || ~ismember(sfs_work(i,6),points_work) || ~ismember(sfs_work(i,10),points_work)
        sfs_work(i,:)=[];
    end
end
sfs_work=sfs_work(:,1);

%突出工作区域
figure(4);
hold on;
for i=ls_notwork'
    plot3([lines0(i,3),lines0(i,7)],[lines0(i,4),lines0(i,8)],[lines0(i,5),lines0(i,9)],'g-');
end
for i=ls_work'
    plot3([lines0(i,3),lines0(i,7)],[lines0(i,4),lines0(i,8)],[lines0(i,5),lines0(i,9)],'b-');
end

clearvars -except vl sfs0 sfs_v0 sfs_work ls_work ls_notwork ps_work lines0 points0 surfaces0 R D d F dk alpha beta r P O_d
