function [ res ] = rmse( X, Y )
    res = sqrt (sum((X-Y).^2)/numel(X));
end

