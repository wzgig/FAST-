%绘制抛物面
x=0;
y=0;
z=-300.0339;
for i=5:5:150
    z=[z;(i^2/559.3357-300.0339)*ones(60,1)];
    for j=0:pi/30:2*pi-0.01
        x=[x;i*cos(j)];
        y=[y;i*sin(j)];
    end
end
dx=[cos(alpha)*cos(beta-pi/2) sin(alpha)*cos(beta-pi/2) sin(beta-pi/2)];
dy=[cos(alpha+pi/2) sin(alpha+pi/2) 0];
dz=[cos(alpha)*cos(beta) sin(alpha)*cos(beta) sin(beta)];
points=[];
for i=1:size(x,1)
    points=[points;x(i)*dx+y(i)*dy+z(i)*dz];
end
plot3(points(:,1),points(:,2),points(:,3),'r-');

sfs0=surfaces0;
sfs_v0=[];%sfs0中，所有平面的方程系数、反射角笛卡尔坐标参数
vl=-P/norm(P);%笛卡尔坐标系入射光源的向量参数
for i=1:size(sfs0,1)
    matrix=[sfs0(i,3) sfs0(i,4) 1;sfs0(i,7) sfs0(i,8) 1;sfs0(i,11) sfs0(i,12) 1];
    solution=matrix\[matrix -[sfs0(i,5);sfs0(i,9);sfs0(i,13)]];%平面方程：ax+by+z+d=0。z方向始终为正，故z方向系数取1
    v_d=[solution(1,end) solution(2,end) 1];%笛卡尔坐标系下法向量的向量参数
    v_d=v_d/norm(v_d);
    v_d=v_d*dot(v_d,vl);
    v_r=[2*v_d(1)-vl(1) 2*v_d(2)-vl(2) 2*v_d(3)-vl(3)];%笛卡尔坐标系下反射光向量参数
    sfs_v0=[sfs_v0;sfs0(i,1) solution(1,end) solution(2,end) 1 solution(3,end) v_r(1) v_r(2) v_r(3)];
end

clearvars -except vl sfs0 sfs_v0 sfs_work ls_work ls_notwork ps_work lines0 points0 surfaces0 R D d F dk alpha beta r P O_d
