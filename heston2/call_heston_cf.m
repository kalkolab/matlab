function y = call_heston_cf(s0, v0, vbar, a, vvol, r, rho, t, k)
% Heston call value using characteristic functions.
% y = call_heston_cf(s0, v0, vbar, a, vvol, r, rho, t, k)
% Inputs:
% s0: stock price
% v0: initial volatility (v0^2 initial variance)
% vbar: long-term variance mean
% a: variance mean-reversion speed
% vvol: volatility of the variance process
% r: risk-free rate
% rho: correlation between the Weiner processes of the stock price and its variance
% t: time to maturity
% k: option strike
% chfun_heston: Heston characteristic function
% 1st step: calculate pi1 and pi2
 % Inner integral 1
int1 = @(w, s0, v0, vbar, a, vvol, r, rho, t, k) real(exp(-i.*w*log(k)).*chfun_heston(s0, v0, vbar, a, vvol, r, ...
rho, t, w-i)./(i*w.*chfun_heston(s0, v0, vbar, a, vvol, r, rho, t, -i))); % inner integral1
int1 = integral(@(w)int1(w,s0, v0, vbar, a, vvol, r, rho, t, k),0,100); % numerical integration
pi1 = int1/pi+0.5; % final pi1
 % Inner integral 2:
int2 = @(w, s0, v0, vbar, a, vvol, r, rho, t, k) real(exp(-i.*w*log(k)).*chfun_heston(s0, v0, vbar, a, vvol, r, ...
rho, t, w)./(i*w));
int2 = integral(@(w)int2(w,s0, v0, vbar, a, vvol, r, rho, t, k),0,100);int2 = real(int2);
pi2 = int2/pi+0.5; % final pi2

% 2rd step: calculate call value
y = s0*pi1-exp(-r*t)*k*pi2;
end
