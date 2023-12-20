clear;
clc;
% 电子回旋频率，单位：MHz
fce = 1.4;
% f2峰值高度，单位：km
hmF2 = 300;
% 电离层下界
HE = 10;
% 电离层上界
Hmax = 700;
% 观测仪器的发射频率，单位：MHz
f = 0.01:0.01:10;

% 构造真高数组，分辨率1000
hrO = zeros(1,1000);
hrX = zeros(1,1000);

% 分别计算O波和X波的临界反射高度
% 假设连续发射1～10MHz的信号
% 当折射指数为0，信号到达反射高度
for m = 1:1000
    % 发射信号的下界和上界 为方便取为fc/2 ~ 等离子体频率最大值（在hmF2处取得）
    if(f(m)>=fce/2.0 && f(m)<=fplasma(hmF2))
        % hpropagate 信号到达的高度，以0.01km的步长向上累积
        hpropagate = 0.01;
        [refrac2O,~,~,~] = refractiveAH(fplasma(hpropagate),fce,f(m));
        % 当折射指数在0到1之间，信号持续传播，高度递增
        while(refrac2O>=0 && hpropagate<=Hmax)
            hpropagate = hpropagate+0.01;
            [refrac2O,~,~,~] = refractiveAH(fplasma(hpropagate),fce,f(m));
        end
        % O波的反射高度，对超出范围进行回退
        hrO(m) = hpropagate-0.01;
    else
        % 发射信号以外的范围，真高设为空值
        hrO(m) = NaN;
    end
    % X波的临界频率略高于O波,因此发射频率控制在fc～fp+fc/2之间
    if(f(m)>=fce && f(m)<=fplasma(hmF2)+fce/2.0+0.02)
        hpropagate = HE;
        [~,refrac2X,~,~] = refractiveAH(fplasma(hpropagate),fce,f(m));
        while(refrac2X>=0 && refrac2X<=1 && hpropagate<=Hmax)
            hpropagate = hpropagate+0.01;
            [~,refrac2X,~,~] = refractiveAH(fplasma(hpropagate),fce,f(m));
        end
        % X波的反射高度，对超出范围进行回退
        hrX(m) = hpropagate-0.01;
    else
        hrX(m) = NaN;
    end
end

% 构造虚高数组，分辨率同样取1000
hvO = ones(1,1000);
hvX = ones(1,1000);

% 利用折射指数对虚高进行积分 原理： hv = integeral(n_group,0,hr) 
% ng 表示群折射指数
for m = 1:1000
    % 真高非空，表示在发射频率范围内
    if(~isnan(hrO(m)))
        % 积分范围 0～反射高度之间
        h = 0:0.01:hrO(m);
        % 计算O波的群折射指数
        [~,~,refracO_group,~] = refractiveAH(fplasma(h),fce,f(m));
        % O波的的传播路径积分，离散值复杂函数用梯形积分 trapz
        hvO(m) = trapz(h,refracO_group);
    else
        % 发射频率范围以外，虚高设为空值
        hvO(m) = NaN;
    end
    % X波的的传播路径积分
    if(~isnan(hrX(m)))
        h = 0:0.01:hrX(m);
        [~,~,~,refracX_group] = refractiveAH(fplasma(h),fce,f(m));
        hvX(m) = trapz(h,refracX_group);
    else
        hvX(m) = NaN;
    end
end

% 绘图
figure(1);
% 电离层有效范围 
h = HE:Hmax;
% 电子等离子体频率剖面
plot(fplasma(h),h,'-',Color=[0.4940 0.1840 0.5560],DisplayName='f_{plasma}');
hold on;
% O波描迹
plot(f,hvO,'.',Color=[0.6350 0.0780 0.1840],DisplayName= 'O 波模');
hold on;
% X波描迹
plot(f,hvX,'.',Color=[0 0.4470 0.7410],DisplayName='X 波模');
% 图例
legend;
% 坐标轴标记
xlabel('f / MHz');
ylabel('h / km');
% 出图
print(gcf,'-dpng','频高图.png');
