function y = chfun_heston(s0, v0, vbar, a, vvol, r, rho, t, w);
% Heston characteristic function.
% Inputs:
% s0: stock price
% v0: initial volatility (v0^2 initial variance)
% vbar: long-term variance mean
% a: variance mean-reversion speed
% vvol: volatility of the variance process
% r : risk-free rate
% rho: correlation between the Weiner processes for the stock price and its variance
% w: points at which to evaluate the function
% Output:
% Characteristic function of log (St) in the Heston model
% Interim calculations
alpha = -w.*w/2 - i*w/2;
beta = a - rho*vvol*i*w;
gamma = vvol*vvol/2;
h = sqrt(beta.*beta - 4*alpha*gamma);
rplus = (beta + h)/vvol/vvol;
rminus = (beta - h)/vvol/vvol;
g=rminus./rplus;
% Required inputs for the characteristic function
C = a * (rminus * t - (2 / vvol^2) .* log((1 - g .* exp(-h*t))./(1-g)));
D = rminus .* (1 - exp(-h * t))./(1 - g .* exp(-h*t));
% Characteristic function evaluated at points w
y = exp(C*vbar + D*v0 + i*w*log(s0*exp(r*t))); 