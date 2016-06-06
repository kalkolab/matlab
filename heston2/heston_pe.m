function [ val flag ] = heston_pe(x, S, r, T, K, M )
    global heston_data;
    flag = 1;
    val = sum((call_heston_v(heston_data(:,1), x(1), x(2),  (x(5)+x(3)^2)/(2*x(2)), x(3), heston_data(:,4), x(4), heston_data(:,3), heston_data(:,2) ) - heston_data(:,5)).^2);
    
end

