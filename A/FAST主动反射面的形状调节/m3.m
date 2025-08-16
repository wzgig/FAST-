isOk=[];
numOk=0;
score=0;
P0=P/norm(P)*R;
for j=1:size(sfs0,1)
    sf_v=sfs_v0(j,2:5);
    v_r=[sfs_v0(j,6:8);P];
    t=(-sf_v(1)*v_r(2,1)-sf_v(2)*v_r(2,2)-sf_v(3)*v_r(2,3)-sf_v(4))/(sf_v(1)*v_r(1,1)+sf_v(2)*v_r(1,2)+sf_v(3)*v_r(1,3));
    p0=v_r(1,:)*t+v_r(2,:);%找到该平面上能反射光线到P点的点
    v1=sfs0(j,3:5)-sfs0(j,7:9);
    v2=sfs0(j,11:13)-sfs0(j,7:9);
    v3=sfs0(j,3:5)-sfs0(j,11:13);
    S0=sin(acos(dot(v1,v2)/(norm(v1)*norm(v2))))*norm(v1)*norm(v2);%求三角形面积，后期用于判断点和三角形的关系
    S=0;
    v=[v1;v2;v3];
    for k=1:3
        m=-v(k,1)*p0(1)-v(k,2)*p0(2)-v(k,3)*p0(3);
        t=(-v(k,1)*sfs0(j,4*k-1)-v(k,2)*sfs0(j,4*k)-v(k,3)*sfs0(j,4*k+1)-m)/(v(k,1)^2+v(k,2)^2+v(k,3)^2);
        p=v(k,:)*t+sfs0(j,4*k-1:4*k+1);
        S=S+norm(p0-p)*norm(v(k,:));
    end
    if S<S0+0.000001
        isOk=[isOk;1];
        numOk=numOk+1;
        score=score+10000;
        score=score-norm(p0-sfs0(j,3:5))-norm(p0-sfs0(j,7:9))-norm(p0-sfs0(j,11:13));
        plot3([P(1),p0(1)],[P(2),p0(2)],[P(3),p0(3)],'k-');
    else
        score=score-100*norm(p0-(sfs0(j,3:5)+sfs0(j,7:9)+sfs0(j,11:13))/3);
        isOk=[isOk;0];
    end
end
scores=[0 score numOk];%历次成绩

clearvars -except numOk isOk scores vl sfs0 sfs_v0 sfs_work ls_work ls_notwork ps_work lines0 points0 surfaces0 R D d F dk alpha beta r P O_d
