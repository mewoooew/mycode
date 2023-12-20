function [refrac2O,refrac2X,refracO_group,refracX_group] = refractiveAH(fpe,fce,femit)
% 根据A-H磁离子理论计算电磁波在等离子体中传播过程中的折射指数

% X = 等离子体频率^2/发射频率^2
X = fpe.^2/femit.^2;
% X = 电子回旋频率^2/发射频率^2
Y = fce/femit;

% O波的折射指数的平方值
refrac2O = 1-X;
% X波的折射指数的平方值
refrac2X = 1-X./(1-Y^2./(1-X));

% O波的群折射指数 nO_group * nO = 1
refracO_group = sqrt(1./refrac2O);
% X波的群折射指数 nX_group * nX ~ 1
refracX_group = sqrt(1./refrac2X);
%refracX_group = (1+X.*Y^2./(X+Y^2-1).^2).*sqrt(1./refrac2X);
end

