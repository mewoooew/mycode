clear;
clc;
% cyctron frequency
fc=1.4;
% 观测仪器的发射频率
f = 0.01:0.01:10;

% O波和X波的A-H关系式（折射指数的平方 = ……），是高度（或者说fp）和发射频率的函数
% ionogram计算对应高度的fp,是高度的函数
% 计算过程中ionogram(x)是数组，因此用到矩阵运算符号
frac2O = @(x,femit) 1-ionogram(x).^2/femit^2;
frac2X = @(x,femit) 1-(ionogram(x).^2/femit^2)./(1-(fc^2/femit^2)./(1-ionogram(x).^2/femit^2));

% 构造虚高数组，分辨率1000
hvO = zeros(1000);
hvX = zeros(1000);

% 分别计算O波和X波的临界反射高度
% 假设连续发射1～10MHz的信号
% 当折射指数在0到1之间，高度一直增加
% 当高度不再增加，则表示信号到达反射高度
for m = 1:1000
    %发射信号的下界和上界 为方便取为fc/2 ~ 等离子体频率最大值（在300km高度处取得）
    if(f(m)>=fc/2&&f(m)<=ionogram(300))
        % k表示信号到达的高度，以0.01km的步长向上累积
        k = 0.01;
        while(frac2O(k,f(m))>0&&k<=600)
            k = k+0.01;
        end
        % O波的临界反射高度
        hvO(m) = k-0.01;
    end
    % X波的临界频率略高于O波,因此发射频率控制在fc～fp+fc/2之间
    if(f(m)>=fc&&f(m)<=ionogram(300)+fc/2)
        k = 0.01;
        while(frac2X(k,f(m))>0&&frac2X(k,f(m))<=1&&k<=600)
            k = k+0.01;
        end
        % X波的临界反射高度
        hvX(m) = k-0.01;
    end
end

% 构造实高数组，分辨率同样取1000
hrO = ones(1000);
hrX = ones(1000);
% 利用折射指数对虚高进行积分 原理： n = c/v
for m = 1:1000
    % 积分范围 0～临界反射高度之间
    h = 0:0.01:hvO(m);
    % O波的的传播路径积分
    hrO(m) = trapz(h,sqrt(1./frac2O(h,f(m))));
    % X波的的传播路径积分
    h = 0:0.01:hvX(m);
    hrX(m) = trapz(h,sqrt(1./frac2X(h,f(m))));
end

% 绘图
figure(1);
h = 50:900;
plot(ionogram(h),h,'-',Color=[0.4940 0.1840 0.5560]);
hold on;
plot(f,hrO,'.',Color=[0.6350 0.0780 0.1840]);
hold on;
plot(f,hrX,'.',Color=[0 0.4470 0.7410]);
legend('f_{plasma}','O描迹','','X描迹');
xlabel('f / MHz');
ylabel('h / km');
f = gcf;
% 出图
print(gcf,'-dpng','频高图.png');
