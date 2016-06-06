function y = call_heston_v(S, v0, vbar, a, vvol, r, rho, T, K )

N = numel(S);


for i=1:N
    y(i,1) = call_heston_cf(S(i), v0, vbar, a, vvol, r(i), rho, T(i), K(i));
end


end

