% Heston Model Calibration
%
% Options2011: set of options with on index
% 
% Author: Jonathan Frei, 2016
%


load('Options2011.mat');
% K is a vector of strikes 
% PC is a vector containing either a 1 or 2 (1 for call, 2 for put)
% Price is a vector of market prices of the options 
% S is a vector of spot prices
% T is a vector of maturities (must be the same unit as r and q)
% dateO is the starting date of the option
% r is a vector of interest rates
% q is a vector of dividend rates



% PC==1 for calls, PC==2 is for puts
% dateO==min(dateO) selects the first dates only (Wednesdays here)
indx = find(dateO==min(dateO)&PC==1); 

% reducing vectors to selected options
K = K(indx);
S = S(indx);
T = T(indx);
t = zeros(size(T)); 
r = r(indx);
q = q(indx);
PC = PC(indx);
MarketPrice = Price(indx);


%MarketIVCheck = blsimpv(S, K, r, T, MarketPrice, [], q, [], (PC==1))*100;
[ MarketIV, dates ] = bsmivec( MarketPrice, S,K,T,t,r, q, PC );



tic
startparameters = [0.2^2 2 -0.6 0.25^2 1  ];
% options = optimoptions('lsqnonlin', 'Display', 'iter');
% xopt = lsqnonlin(@(x) Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice,...
%     startparameters, ... %start values
%     ... % v0,theta,rho,kappa,sigma
%     [eps eps -1+eps eps eps  ], ... % lower bound for parameter vector
%     [Inf Inf 1-eps Inf  Inf  ], ...  % upper bound for parameter vector
%     options);

%[xopt,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,startparameters,[eps eps -1+eps eps eps  ],[Inf Inf 1-eps Inf  Inf  ])
ObjectivePE = @(x) sum((Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice).^2);
ObjectiveVolErr = @(x) sum((bsmivec( Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,x(1),x(2),x(3),x(4),x(5)), ...
     S,K,T,t,r, q, PC) - MarketIV).^2);
 ObjectiveAPE = @(x) sum(abs(Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice));
 ObjectiveRPE = @(x) sum(abs(Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice)./MarketPrice);
% 
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

%[xopt,fval,exitFlag,output] = fminunc(ObjectiveAPE, ...
%    startparameters, options)

% 
% % options = optimoptions('particleswarm', 'SwarmSize',50,'HybridFcn',@fmincon, 'Display', 'iter', 'UseParallel', true);
%[xopt,fval,exitFlag,output] = fminsearch(ObjectiveAPE, ...
%    startparameters)

% 
% options = optimoptions('particleswarm', 'SwarmSize',50,'HybridFcn',@patternsearch, 'Display', 'iter', 'UseParallel', true);
% [xopt,fval,exitFlag,output] = particleswarm(ObjectivePE, ...
%     5, ...
%     [eps eps -1+eps eps eps  ],...
%     [Inf Inf 1-eps Inf  Inf  ], ...
%     options)
options = saoptimset('PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
[xopt,fval,exitFlag,output] = simulannealbnd(ObjectivePE, ...
    startparameters, ...
    [eps eps -1+eps eps eps  ],...
    [Inf Inf 1-eps Inf  Inf], ...
    options);

% options = optimoptions('particleswarm', 'SwarmSize',50,'HybridFcn',@fmincon, 'Display', 'iter', 'UseParallel', true);
% [xopt,fval,exitFlag,output] = ga(ObjectiveFunction, ...
%     5)


toc



disp(['Optimal parameter vector: ' num2str(xopt)]);
ModelPrices =  Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,xopt(1),xopt(2),xopt(3),xopt(4),xopt(5))
%ModelIVCheck = blsimpv(S, K, r, T, ModelPrices, [], q, [], (PC==1));
ModelIV = bsmivec( ModelPrices, S,K,T,t,r, q, PC );
disp(['RMSE: ' num2str(sqrt(mean((ModelPrices-MarketPrice).^2)))]);
figure;
DifferentT = sort(unique(T));
for i=1:numel(DifferentT)
   subplot(ceil(numel(unique(T))/2),2,i);
   indx = (T==DifferentT(i));
   %plot(K(indx)./S(indx), MarketIVCheck(indx), 'Marker', 'o', 'LineStyle', 'no');
   plot(K(indx)./S(indx), MarketIV(indx), 'Marker', 'o', 'LineStyle', 'no');
   title([num2str(round(DifferentT(i)*365)) ' days to maturity']);
   hold on; 
   plot(K(indx)./S(indx), ModelIV(indx), 'Marker', 'x', 'LineStyle', 'no', 'Color', 'r');
   %plot(K(indx)./S(indx), ModelIVCheck(indx), 'Marker', 'x', 'LineStyle', 'no', 'Color', 'r');
   legend({'Market', 'Model'});
   xlabel('Moneyness K/S'); ylabel('implied volatility');
end


