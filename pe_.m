function [ error ] = pe_(model, mkt)
%RMSE Summary of this function goes here
%   Detailed explanation goes here
    
N = length(model);
error = 0;
for i=1:N
    error = error + (model(i) - mkt(i))^2;
end
error = sqrt(error/N);
end

