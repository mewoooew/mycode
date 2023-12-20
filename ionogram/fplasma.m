function fplasma = fplasma(height)
% 计算等离子体频率的函数，自变量为高度
% 给定参数
hmF2 = 300;
HF = 50;
hmE = 120;
HE = 10;

% 电子数密度表达式
Ne = 1e12*exp(1-(height-hmF2)/HF-exp((hmF2-height)/HF))+5e10*exp(1-(height-hmE)/HE-exp((hmE-height)/HE));
% 电子数密度到等离子体频率的换算
fplasma = (8.98e-6)*sqrt(Ne);

end
