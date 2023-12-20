clc;
clear;
syms ionogram fc femit;
ftest(ionogram,fc,femit) = 1-(ionogram^2/femit^2-ionogram^4/femit^4)/(1-ionogram^2/femit^2-fc^2/femit^2);
diff(ftest(ionogram,fc,femit),femit);
fnew(ionogram,fc,femit) = ftest(ionogram,fc,femit)+femit/2*diff(ftest(ionogram,fc,femit),femit)
