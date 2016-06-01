syms h T k r w alpha beta gamma rho
a = .5*(k^2-k);
b = beta + alpha*(k-gamma)^2;
B0 = a*(1-b^T)/(1-b);
A0 = T*k*r + (w+alpha)/(1-b)*(a*T-B0);

C0 = A0 + h*B0;

kappa1 = subs(subs(diff(C0, k), k, 0), alpha*gamma^2 + beta, rho)
kappa3 = subs(subs(diff(C0, k, 3), k, 0), alpha*gamma^2 + beta, rho)
kappa4 = subs(subs(diff(C0, k, 4), k, 0), alpha*gamma^2 + beta, rho)