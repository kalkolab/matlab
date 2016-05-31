function [ error ] = pe(price, c, S, X, T, rf, vol, model )
%RMSE Summary of this function goes here
%   Detailed explanation goes here
    
N = length(S);
error = 0;
for i=1:N
    error = error + (price(S(i),X(i),T(i),rf(i),vol(i),model) - c(i))^2;
end
error = sqrt(error/N);
end

