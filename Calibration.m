data = dlmread('sample.csv', ';', 1, 0);
N = size(data,1);

mkt_call = data(:,1);
S = data(:,2);
X = data(:,3);
T = data(:,4);
vol = data(:,5);
r = data(:,6);


%ObjectiveFunction = @(alpha, beta, gamma, lambda, omega) pe(@HNCall2, mkt_call, S,X,T,vol,r, struct('alpha', alpha, 'beta', beta,'gamma', gamma,'lambda', lambda,'omega', omega));
ObjectiveFunction = @(x) pe(@HNCall2, mkt_call, S,X,T,vol,r, struct('alpha', x(1), 'beta', x(2),'gamma', x(3),'lambda', x(4),'omega', x(5)));
%ObjectiveFunction(0.00004,.159, 430, 196,1e-5) 
X0 = [0.00004 .159 430 196 1e-5];
%[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0)
[x,fval,exitflag] = fminsearch(ObjectiveFunction,X0, optimset('PlotFcns',@optimplotx))