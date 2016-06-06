%% load
clear all
global heston_data;

data = dlmread('sample.csv', ';', 1, 0);

% PC==1 for calls, PC==2 is for puts
% dateO==min(dateO) selects the first dates only (Wednesdays here)

% reducing vectors to selected options
data = data(data(:,4) < 365 & data(:,4)>30,:);
K = data(:,3);
S = data(:,2);
T = data(:,4)./365;
t = zeros(size(T)); 
r = ones(size(T)).*0.005;
q = zeros(size(r));
PC = ones(size(S));
MarketPrice = data(:,1);

%call_heston_v(S,K,T-t, r, 0.04, 0.2, 2e-5, 0.01, -0.7)

%% test
%mazzoniCall(S,K,T,t,r-q, 1e-6,.5e-5, 0.6, 100, 0.04)
%% marketIV
MarketIVCheck = blsimpv(S, K, r, T, MarketPrice, [], q, [], (PC==1))*100;
[ MarketIV, dates ] = bsmivec( MarketPrice, S,K,T,t,r, q, PC );

%% calibrate
x0 = [0.16 0.16 1 -0.52 5];
lb = [eps eps eps -1+eps eps ]; % lower bound for parameter vector 
ub = [1 1 6 1 20];

%% lsq
tic


options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunEvals', Inf, 'MaxIter', Inf, 'TolFun', 1e-6);
xopt = lsqnonlin(@(x) call_heston_v(S, x(1), x(2),  (x(5)+x(3)^2)/(2*x(2)), x(3), r, x(4), T-t, K ) - MarketPrice,...
    x0, ... %start values
    ... % v0,theta,rho,kappa,sigma
    lb, ... % lower bound for parameter vector
    ub, ...  % upper bound for parameter vector
    options);
    
% options = optimoptions('fminsearch', 'Display', 'iter', 'MaxFunEvals', Inf, 'MaxIter', Inf, 'TolFun', 1e-6);
% xopt = fminsearch(@(x) sum((mazzoniCall(S,K,T,t,r-q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice).^2),...
%     startparameters, ... %start values
%     options);

elapsedLSQ = toc

%% asamin
heston_data(:,1)=S;
heston_data(:,2)=K;
heston_data(:,3)=T-t;
heston_data(:,4)=r;
heston_data(:,5)=MarketPrice;

tic
asamin('set','rand_seed',696969);
asamin('set','asa_out_file','asatest1.log');
asamin('set','test_in_cost_func',0);
[fstar, xstar, grad, hessian, state] = ...
  asamin('minimize', 'heston_pe', x0', lb', ub', -1*ones(size(x0'))); 

elapsedtime = toc

%% show res
disp(['Optimal parameter vector: ' num2str(xopt)]);
ModelPrices =  mazzoniCall(S,K,T,t,r-q,xopt(1),xopt(2),xopt(3),xopt(4),xopt(5))
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


%% RMSE
ModelPrices =  call_heston_v(S, xstar(1), xstar(2),  (xstar(5)+xstar(3)^2)/(2*xstar(2)), xstar(3), r, xstar(4), T-t, K );
ModelPrices2 =  call_heston_v(S, xopt(1), xopt(2),  (xopt(5)+xopt(3)^2)/(2*xopt(2)), xopt(3), r, xopt(4), T-t, K );

rmse(ModelPrices, MarketPrice)
rmse(ModelPrices2, MarketPrice)