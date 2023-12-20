%% 
% 清空
clc;
clear all;
close all;
tic;  %计时
%%

%%
% 积分范围初始化
x = -4:0.08:4;   %x轴 分辨率
z = -4:0.08:4;   %z轴
xpoints = size(x,2);
zpoints = size(z,2);
% 水平方位角phi的分点
ppoints = 500;
phi = linspace(0,pi,ppoints);   % phi的范围（0,pi）
phi = phi';   %phi向量化，方便对phi积分

% 磁场分量的微分
% 初始化空矩阵
d_Bx = zeros(zpoints,xpoints,ppoints);
d_Bz = zeros(zpoints,xpoints,ppoints);
for P = 1:ppoints
    for n = 1:zpoints
        for m = 1:xpoints
            R = 1+z(n)^2+x(m)^2-2*x(m)*cos(phi(P));
            d_Bx(n,m,P) = z(n)/pi*cos(phi(P))/(R^(3/2));
            d_Bz(n,m,P) = 1.0/pi*(1-x(m)*cos(phi(P)))/(R^(3/2));
        end
    end
end
%%

%%
% 进行积分求出磁场x、z分量
Bx = zeros(zpoints,xpoints);
Bz = zeros(zpoints,xpoints);
for n = 1:zpoints
    for m = 1:xpoints
        Bx(n,m) = trapz(phi,d_Bx(n,m,:));    % 梯形法积分x分量
        Bz(n,m) = trapz(phi,d_Bz(n,m,:));    % z分量
    end
end
alpha = atan2d(Bz,Bx);      % 磁场方向
Bf = sqrt(Bx.^2+Bz.^2);     % 总磁场强度
toc;
%%

%% 
% 磁场分布图
[X,Z] = meshgrid(x,z);    % 建立坐标矩阵

p1 = figure(1);

% z分量图
subplot(2,2,1);
s1 = surf(X,Z,Bz);
zlim([-2,2]);  % 坐标轴范围
pos = axis;    % 坐标轴位置
xlabel("x / a",'Position',[(pos(1)+pos(2))/2.0 pos(3)-2 pos(5)]);   % 标签
ylabel("z / a",'Position',[pos(1)-2 (pos(3)+pos(4))/2.0 pos(5)]);
zlabel("B_{z} / B_{0}");
title("\fontsize{16} 图 1 环电流磁感应强度 z 分量分布面",'Position',[pos(1)-2.5 pos(3)-4 pos(5)]);
set(gca,'ZTick',-2:0.5:2);   % 坐标轴标注间隔
colormap(gca,'cool')
s1.FaceAlpha = 0.9;          % 透明度
s1.EdgeColor = 'none';       % 不绘制边缘
s1.FaceColor = 'interp';     % 绘图上色进行插值平滑
hold on;
% 在图上绘制轴线
equi = NaN(size(Bz));       % 轴线上的变量
equi((zpoints+1)/2,:,:) = Bz((zpoints+1)/2,:,:);
equi(:,(xpoints+1)/2,:) = Bz(:,(xpoints+1)/2,:);
surf(X,Z,equi,EdgeColor="black",LineWidth=2);

% x分量图
subplot(2,2,2);
s2 = surf(X,Z,Bx);
zlim([-2,2]);
pos = axis;
xlabel("x / a",'Position',[(pos(1)+pos(2))/2.0 pos(3)-2 pos(5)]);
ylabel("z / a",'Position',[pos(1)-2 (pos(3)+pos(4))/2.0 pos(5)]);
zlabel("B_{x} / B_{0}");
title("\fontsize{16} 图 2 环电流磁感应强度 x 分量分布面",'Position',[pos(1)-2.5 pos(3)-4 pos(5)]);
set(gca,'ZTick',-2:0.5:2);
colormap(gca,'cool');
s2.FaceAlpha = 0.9;
s2.EdgeColor = 'none';
s2.FaceColor = 'interp';
hold on;
% 在图上绘制轴线
equi = NaN(size(Bx));
equi((zpoints+1)/2,:,:) = Bx((zpoints+1)/2,:,:);
equi(:,(xpoints+1)/2,:) = Bx(:,(xpoints+1)/2,:);
surf(X,Z,equi,EdgeColor="black",LineWidth=2);

% 总强
subplot(2,2,3);
s3 = surf(X,Z,Bf);
zlim([0,2]);
pos = axis;
xlabel("x / a",'Position',[(pos(1)+pos(2))/2.0 pos(3)-2 pos(5)]);
ylabel("z / a",'Position',[pos(1)-2 (pos(3)+pos(4))/2.0 pos(5)]);
zlabel("B / B_{0}");
title("\fontsize{16} 图 3 环电流合磁感应强度分布面",'Position',[pos(1)-2.5 pos(3)-4 pos(5)]);
set(gca,'ZTick',0:0.5:2);
colormap(gca,'winter');
s3.FaceAlpha = 0.85;
s3.EdgeColor = 'none';
s3.FaceColor = 'interp';
% 在图上绘制轴线
hold on;
equi = NaN(size(Bf));
equi((zpoints+1)/2,:,:) = Bf((zpoints+1)/2,:,:);
equi(:,(xpoints+1)/2,:) = Bf(:,(xpoints+1)/2,:);
surf(X,Z,equi,EdgeColor="black",LineWidth=2);

% 磁场方向
subplot(2,2,4);
s4 = surf(X,Z,alpha);
zlim([-200,200]);
pos = axis;
xlabel("x / a",'Position',[(pos(1)+pos(2))/2.0 pos(3)-2 pos(5)]);
ylabel("z / a",'Position',[pos(1)-2 (pos(3)+pos(4))/2.0 pos(5)]);
zlabel("\alpha (\circ)");
title("\fontsize{16} 图 4 环电流磁感应强度方向分布面",'Position',[pos(1)-2.3 pos(3)-2.3 pos(5)]);
set(gca,'ZTick',-200:100:200);
colormap(gca,'cool');
s4.FaceAlpha = 0.8;
s4.EdgeColor = 'none';
s4.FaceColor = 'interp';
% 在图上绘制轴线
hold on;
equi = NaN(size(alpha));
equi((zpoints+1)/2,:,:) = alpha((zpoints+1)/2,:,:);
equi(:,(xpoints+1)/2,:) = alpha(:,(xpoints+1)/2,:);
surf(X,Z,equi,EdgeColor="black",LineWidth=2);
view([-45,60])

%print(p1,'-dpng','-r250','B_field_RingCurrent') %打印
%%

%% 
% 磁场流线图
p2 = figure(2);
[startx,startz] = meshgrid(0:0.1:0.4,0);
h1 = streamline(X,Z,Bx,Bz,startx,startz,[0.2,2000]);
% 坐标轴设置
pos = axis;
xlabel("x / a");%,'Position',[(pos(1)+pos(2))/2.0 pos(3)-2 pos(5)]);
ylabel("z / a");%,'Position',[pos(1)-2 (pos(3)+pos(4))/2.0 pos(5)]);
title("\fontsize{16} 图 5 环电流磁感应线",'Position',[pos(1) pos(3)-3.6]);
% 流线参数设置
lennn = length(h1);
for kkk = 1:lennn
    h1(kkk).Color = [0 kkk/(2*lennn-1)/1.5 kkk/(2*lennn-1)];
    h1(kkk).LineWidth = 2;
end
% 对图像沿y轴和x轴进行翻转复制
h2 = copyobj(h1,gca);
rotate(h2,[1,0,0],180,[0 0 0]);
h3 = copyobj(allchild(gca),gca);
rotate(h3,[0,1,0],180,[0 0 0]);
for kk = 1:(lennn-1)
    [startx,startz] = meshgrid(0.4+kk*0.1,0);
    h1 = streamline(X,Z,Bx,Bz,startx,startz,[0.05/(kk+1),4500]);
    h2 = streamline(-X,Z,-Bx,Bz,-startx,startz,[0.05/(kk+1),4500]);
    h1.Color = [0 (lennn + kk)/(2*lennn-1)/1.5 (lennn + kk)/(2*lennn-1)];
    h1.LineWidth = 2;
    h2.Color = h1.Color;
    h2.LineWidth = 2;
end
% 绘制电流源
radi = 0.1;
rectangle('Position',[1-radi,0-radi,2*radi,2*radi],'Curvature',[1,1]);
rectangle('Position',[-1-radi,0-radi,2*radi,2*radi],'Curvature',[1,1]);
hold on;

% 绘制电流方向
scatter(1,0,180,'black',"X","LineWidth",4);
hold on;
scatter(-1,0,180,'black',".","LineWidth",4);
hold on;
[lrope,wrope] = meshgrid(-1+radi*1.2:0.01:1-radi*1.2,0);
scatter(lrope,wrope,'black',"_","LineWidth",40);
xlim([-4 4]);
ylim([-3 3]);
box on;

%print(p2,'-dpng','-r250','B_stream_RingCurrent')  %打印
%%