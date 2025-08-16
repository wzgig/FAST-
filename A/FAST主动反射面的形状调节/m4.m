%优化计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
times=100;              %优化次数
score_lr=10000000;   %两主索节点间距离判分系数
score_ok=10000;          %期望工作口径系数
score_l=100;       %期望工作口径判分系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迭代的数据：ps0,ls0,sfs0,sfs_v0,par
%参与运算的数据：ps,ls,sfs,sfs_v

%数据介绍（按列）：
%ps0：1_节点序数 2_x坐标 3_y坐标 4_z坐标 5~10_该节点构成的线的序数 11~16_该节点构成的面的序数
%ls0：1_线序数 2_一个端部节点序数 3_该节点x坐标 4_该节点y坐标 5_该节点z坐标 6~9_另一节点信息 10_线段长度
%sfs0：1_面序数 2_一个角节点序数 3_该节点x坐标 4_该节点y坐标 5_该节点z坐标 6~13_另两个节点信息
%sfs_v0：1_面序数 2~5_面的参数方程系数abcd(ax+by+cz+d=0) 6~8_反射光方向单位向量
%par：1_每个节点下一次变化量的大小(伸长或缩短) 2_已变化量

ps0=points0;
ps=ps0;
ls0=lines0;
ls=ls0;
sfs0=surfaces0;
sfs=sfs0;
sfs_v=sfs_v0;
par=[ones(size(ps,1),1)*0.001048576 zeros(size(ps,1),1)];%每次伸缩长度、下拉索已伸缩长度
for time=1:times
    time%实时显示迭代次数
    for i=1:size(ps,1)
        for n=1:2
            l=norm(ps(i,2:4));
            if n==1
                ps(i,2:4)=ps(i,2:4)/l*(l-par(ps(i,1),1));%先减后加运算速度加快
            elseif n==2
                ps(i,2:4)=ps(i,2:4)/l*(l+par(ps(i,1),1));
            end
            if abs(norm(ps(i,2:4))-norm(points0(i,2:4)))<=0.6
                sco0=0;
                sco=0;
                P0=P/norm(P)*R;
                for j=5:10
                    if ps(i,j)~=0
                        if ls(ps(i,j),2)==ps(i,1)
                            ls(ps(i,j),3:5)=ps(i,2:4);
                            ls(ps(i,j),10)=norm(ls(ps(i,j),7:9)-ls(ps(i,j),3:5));
                            l_d0=abs(ls0(ps(i,j),10)-lines0(ps(i,j),10))/lines0(ps(i,j),10);
                            l_d=abs(ls(ps(i,j),10)-lines0(ps(i,j),10))/lines0(ps(i,j),10);
                            if l_d0>0.0007
                                p=(ls0(ps(i,j),3:5)+ls0(ps(i,j),7:9))/2;
                                sco0=sco0-(l_d-0.00069999)*score_lr*norm(p-P0)*time;
                            end
                            if l_d>0.0007
                                p=(ls(ps(i,j),3:5)+ls(ps(i,j),7:9))/2;
                                sco=sco-(l_d-0.00069999)*score_lr*norm(p-P0)*time;
                            end
                        elseif ls(ps(i,j),5)==ps(i,1)
                            ls(ps(i,j),7:9)=ps(i,2:4);
                            ls(ps(i,j),10)=norm(ls(ps(i,j),7:9)-ls(ps(i,j),3:5));
                            l_d0=abs(ls0(ps(i,j),10)-lines0(ps(i,j),10))/lines0(ps(i,j),10);
                            l_d=abs(ls(ps(i,j),10)-lines0(ps(i,j),10))/lines0(ps(i,j),10);
                            if l_d>0.0007
                                p=(ls0(ps(i,j),3:5)+ls0(ps(i,j),7:9))/2;
                                sco0=sco0-(l_d-0.00069999)*score_lr*norm(p-P0)*time;
                            end
                            if l_d>0.0007
                                p=(ls(ps(i,j),3:5)+ls(ps(i,j),7:9))/2;
                                sco=sco-(l_d-0.00069999)*score_lr*norm(p-P0)*time;
                            end
                        end
                    end
                end
                for j=11:16
                    if ps(i,j)~=0
                        if sfs(ps(i,j),2)==ps(i,1)
                            sfs(ps(i,j),3:5)=ps(i,2:4);
                            matrix=[sfs(ps(i,j),3) sfs(ps(i,j),4) 1;sfs(ps(i,j),7) sfs(ps(i,j),8) 1;sfs(ps(i,j),11) sfs(ps(i,j),12) 1];
                            solution=matrix\[matrix -[sfs(ps(i,j),5);sfs(ps(i,j),9);sfs(ps(i,j),13)]];%平面方程：ax+by+z+d=0。z方向始终为正，故z方向系数取1
                            v_d=[solution(1,end) solution(2,end) 1];%笛卡尔坐标系下法向量的向量参数
                            v_d=v_d/norm(v_d);
                            v_d=v_d*dot(v_d,vl);
                            v_r=[2*v_d(1)-vl(1) 2*v_d(2)-vl(2) 2*v_d(3)-vl(3)];%笛卡尔坐标系下反射光向量参数
                            sfs_v(ps(i,j),:)=[ps(i,j) solution(1,end) solution(2,end) 1 solution(3,end) v_r(1) v_r(2) v_r(3)];
                        elseif sfs(ps(i,j),6)==ps(i,1)
                            sfs(ps(i,j),7:9)=ps(i,2:4);
                            matrix=[sfs(ps(i,j),3) sfs(ps(i,j),4) 1;sfs(ps(i,j),7) sfs(ps(i,j),8) 1;sfs(ps(i,j),11) sfs(ps(i,j),12) 1];
                            solution=matrix\[matrix -[sfs(ps(i,j),5);sfs(ps(i,j),9);sfs(ps(i,j),13)]];%平面方程：ax+by+z+d=0。z方向始终为正，故z方向系数取1
                            v_d=[solution(1,end) solution(2,end) 1];%笛卡尔坐标系下法向量的向量参数
                            v_d=v_d/norm(v_d);
                            v_d=v_d*dot(v_d,vl);
                            v_r=[2*v_d(1)-vl(1) 2*v_d(2)-vl(2) 2*v_d(3)-vl(3)];%笛卡尔坐标系下反射光向量参数
                            sfs_v(ps(i,j),:)=[ps(i,j) solution(1,end) solution(2,end) 1 solution(3,end) v_r(1) v_r(2) v_r(3)];
                        elseif sfs(ps(i,j),10)==ps(i,1)
                            sfs(ps(i,j),11:13)=ps(i,2:4);
                            matrix=[sfs(ps(i,j),3) sfs(ps(i,j),4) 1;sfs(ps(i,j),7) sfs(ps(i,j),8) 1;sfs(ps(i,j),11) sfs(ps(i,j),12) 1];
                            solution=matrix\[matrix -[sfs(ps(i,j),5);sfs(ps(i,j),9);sfs(ps(i,j),13)]];%平面方程：ax+by+z+d=0。z方向始终为正，故z方向系数取1
                            v_d=[solution(1,end) solution(2,end) 1];%笛卡尔坐标系下法向量的向量参数
                            v_d=v_d/norm(v_d);
                            v_d=v_d*dot(v_d,vl);
                            v_r=[2*v_d(1)-vl(1) 2*v_d(2)-vl(2) 2*v_d(3)-vl(3)];%笛卡尔坐标系下反射光向量参数
                            sfs_v(ps(i,j),:)=[ps(i,j) solution(1,end) solution(2,end) 1 solution(3,end) v_r(1) v_r(2) v_r(3)];
                        end
                        %%%%%%%%%%%%%%%%%%%%
                        %算原来的分数
                        sf_v0=sfs_v0(ps(i,j),2:5);
                        v_r0=[sfs_v0(ps(i,j),6:8);P];
                        t0=(-sf_v0(1)*v_r0(2,1)-sf_v0(2)*v_r0(2,2)-sf_v0(3)*v_r0(2,3)-sf_v0(4))/(sf_v0(1)*v_r0(1,1)+sf_v0(2)*v_r0(1,2)+sf_v0(3)*v_r0(1,3));
                        p0=v_r0(1,:)*t0+v_r0(2,:);%找到该平面上能反射光线到P点的点
                        v1=sfs0(ps(i,j),3:5)-sfs0(ps(i,j),7:9);
                        v2=sfs0(ps(i,j),11:13)-sfs0(ps(i,j),7:9);
                        v3=sfs0(ps(i,j),3:5)-sfs0(ps(i,j),11:13);
                        S0=sin(acos(dot(v1,v2)/(norm(v1)*norm(v2))))*norm(v1)*norm(v2);%求三角形面积，后期用于判断点和三角形的关系
                        S=0;
                        v=[v1;v2;v3];
                        for k=1:3
                            m=-v(k,1)*p0(1)-v(k,2)*p0(2)-v(k,3)*p0(3);
                            t0=(-v(k,1)*sfs0(ps(i,j),4*k-1)-v(k,2)*sfs0(ps(i,j),4*k)-v(k,3)*sfs0(ps(i,j),4*k+1)-m)/(v(k,1)^2+v(k,2)^2+v(k,3)^2);
                            p=v(k,:)*t0+sfs0(ps(i,j),4*k-1:4*k+1);
                            S=S+norm(p0-p)*norm(v(k,:));
                        end
                        if S<S0+0.00000001
                            sco0=sco0+score_ok;
                            sco0=sco0-norm(p0-(sfs0(ps(i,j),3:5)+sfs0(ps(i,j),7:9)+sfs0(ps(i,j),11:13))/3)^2;
                        else
                            sco0=sco0-score_l*norm(p0-(sfs0(ps(i,j),3:5)+sfs0(ps(i,j),7:9)+sfs0(ps(i,j),11:13))/3);
                        end
                        %%%%%%%%%%%%%%%%%%%%
                        %算改变后的分数
                        sf_v=sfs_v(ps(i,j),2:5);
                        v_r=[sfs_v(ps(i,j),6:8);P];
                        t=(-sf_v(1)*v_r(2,1)-sf_v(2)*v_r(2,2)-sf_v(3)*v_r(2,3)-sf_v(4))/(sf_v(1)*v_r(1,1)+sf_v(2)*v_r(1,2)+sf_v(3)*v_r(1,3));
                        p0=v_r(1,:)*t+v_r(2,:);%找到该平面上能反射光线到P点的点
                        v1=sfs(ps(i,j),3:5)-sfs(ps(i,j),7:9);
                        v2=sfs(ps(i,j),11:13)-sfs(ps(i,j),7:9);
                        v3=sfs(ps(i,j),3:5)-sfs(ps(i,j),11:13);
                        S0=sin(acos(dot(v1,v2)/(norm(v1)*norm(v2))))*norm(v1)*norm(v2);%求三角形面积，后期用于判断点和三角形的关系
                        S=0;
                        v=[v1;v2;v3];
                        for k=1:3
                            m=-v(k,1)*p0(1)-v(k,2)*p0(2)-v(k,3)*p0(3);
                            t=(-v(k,1)*sfs(ps(i,j),4*k-1)-v(k,2)*sfs(ps(i,j),4*k)-v(k,3)*sfs(ps(i,j),4*k+1)-m)/(v(k,1)^2+v(k,2)^2+v(k,3)^2);
                            p=v(k,:)*t+sfs(ps(i,j),4*k-1:4*k+1);
                            S=S+norm(p0-p)*norm(v(k,:));
                        end
                        if S<S0+0.00000001
                            sco=sco+score_ok;
                            sco=sco-norm(p0-(sfs(ps(i,j),3:5)+sfs(ps(i,j),7:9)+sfs(ps(i,j),11:13))/3)^2;
                        else
                            sco=sco-score_l*norm(p0-(sfs(ps(i,j),3:5)+sfs(ps(i,j),7:9)+sfs(ps(i,j),11:13))/3);
                        end
                        %%%%%%%%%%%%%%%%%%%%
                    end
                end
                %迭代的数据：ps0,ls0,sfs0,sfs_v0,par
                %参与运算的数据：ps,ls,sfs,sfs_v
                if sco-sco0<0
                    ps=ps0;
                    ls=ls0;
                    sfs=sfs0;
                    sfs_v=sfs_v0;
                else
                    ps0=ps;
                    ls0=ls;
                    sfs0=sfs;
                    sfs_v0=sfs_v;
                    par(i,2)=norm(ps0(i,2:4))-norm(points0(i,2:4));
                    break;
                end
            else
                ps=ps0;
            end
            if n==2 && par(ps(i,1),1)>0.00000001
                par(ps(i,1),1)=par(ps(i,1),1)/2;
            end
        end
    end
    isOk=[];
    numOk=0;
    score=0;
    if time/100-fix(time/100)==0
        figure(time/100+5);
        hold on;
        for i=ls_notwork'
            plot3([ls0(i,3),ls0(i,7)],[ls0(i,4),ls0(i,8)],[ls0(i,5),ls0(i,9)],'g-');
        end
        for i=ls_work'
            plot3([ls0(i,3),ls0(i,7)],[ls0(i,4),ls0(i,8)],[ls0(i,5),ls0(i,9)],'b-');
        end
    end
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
        if S<S0+0.00000001
            isOk=[isOk;1];
            numOk=numOk+1;
            score=score+score_ok;
            score=score-norm(p0-sfs0(j,3:5))-norm(p0-sfs0(j,7:9))-norm(p0-sfs0(j,11:13));
            if time/100-fix(time/100)==0
                plot3([P(1),p0(1)],[P(2),p0(2)],[P(3),p0(3)],'k-');
            end
        else
            score=score-score_l*norm(p0-(sfs0(j,3:5)+sfs0(j,7:9)+sfs0(j,11:13))/3);
            isOk=[isOk;0];
        end
    end
    scores=[scores;time score numOk];%历次成绩
    if time/100-fix(time/100)==0
        R=300;
    end
end