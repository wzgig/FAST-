clc;
clear;
load rawData;

%点坐标连线
figure(1);
plot3(points0(:,2),points0(:,3),points0(:,4));

%FAST全貌

figure(2);
hold on;
for i=1:size(lines0,1)
    plot3([lines0(i,3),lines0(i,7)],[lines0(i,4),lines0(i,8)],[lines0(i,5),lines0(i,9)],'g-');
end

%观测每两根线之间的距离
figure(3);
plot(1:size(lines0,1),lines0(:,10));

clearvars -except lines0 points0 surfaces0
