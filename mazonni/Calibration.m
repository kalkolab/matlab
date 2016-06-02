%% load
data = dlmread('sample.csv', ';', 1, 0);

% PC==1 for calls, PC==2 is for puts
% dateO==min(dateO) selects the first dates only (Wednesdays here)

% reducing vectors to selected options
K = data(:,3);
S = data(:,2);
T = data(:,4);
t = zeros(size(T)); 
r = data(:,6);
q = zeros(size(r));
PC = ones(size(S));
MarketPrice = data(:,1);


%% test
%mazzoniCall(S,K,T,t,r-q, 1e-6,.5e-5, 0.6, 100, 0.04)
%% marketIV
MarketIVCheck = blsimpv(S, K, r, T, MarketPrice, [], q, [], (PC==1))*100;
[ MarketIV, dates ] = bsmivec( MarketPrice, S,K,T,t,r, q, PC );

%% calibrate

tic

startparameters = [1e-6 .5e-5 0.6 100 0.04 ];
options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunEvals', Inf, 'MaxIter', Inf, 'TolFun', 1e-6);
xopt = lsqnonlin(@(x) mazzoniCall(S,K,T,t,r-q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice,...
    startparameters, ... %start values
    ... % v0,theta,rho,kappa,sigma
    [eps eps eps eps eps ], ... % lower bound for parameter vector
    [Inf Inf Inf  Inf Inf ], ...  % upper bound for parameter vector
    options);

% options = optimoptions('fminsearch', 'Display', 'iter', 'MaxFunEvals', Inf, 'MaxIter', Inf, 'TolFun', 1e-6);
% xopt = fminsearch(@(x) sum((mazzoniCall(S,K,T,t,r-q,x(1),x(2),x(3),x(4),x(5)) - MarketPrice).^2),...
%     startparameters, ... %start values
%     options);

toc


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


