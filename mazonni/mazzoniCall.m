function [ prices ] = mazzoniCall( S, K, T, t, mu, w, alpha, beta, gamma, sigma, r )
%MAZONNI Summary of this function goes here
%   Detailed explanation goes here
    % force column vector
    S=S(:);
    K=K(:);
    T=T(:);
    t=t(:);
    mu=mu(:);

    nos = numel(S);
    prices=NaN(nos,1);
    tau=T-t;
    
    F = S.*exp(mu.*tau);

    h = var(S);
    
    rho = alpha*gamma^2 + beta;
    
    d2 = (log(S./K) + mu)./sigma;
    d1 = d2 + sigma;
    
    %c3 = (3*alpha*h*(rho.^T - 1))/(rho - 1)^2 - (6*alpha*(T/2 - (rho.^T - 1)/(2*(rho - 1)))*(alpha + w))/(rho - 1)^2 - ((alpha + w)*((12*alpha^2*gamma^2*(rho.^T - 1))/(rho - 1)^3 - (3*alpha*(rho.^T - 1))/(rho - 1)^2 + (3*T.*alpha.*rho.^(T - 1))/(rho - 1) - (6*alpha*gamma*(rho.^T - 1))./(rho - 1)^2 + (6*T.*alpha*gamma.*rho.^(T - 1))./(rho - 1) - (12*T.*alpha^2*gamma^2*rho.^(T - 1))./(rho - 1)^2 + (6*T.*alpha^2*gamma^2.*rho.^(T - 2).*(T - 1))./(rho - 1)))./(rho - 1) - (12*alpha^2*gamma^2*h*(rho.^T - 1))./(rho - 1)^3 - (3*T.*alpha*h.*rho.^(T - 1))./(rho - 1) + (6*alpha*gamma*h*(rho.^T - 1))./(rho - 1)^2 - (6*alpha*gamma*(alpha + w)*(T - (rho.^T - 1)./(rho - 1) + (2*alpha*gamma*(rho.^T - 1))./(rho - 1)^2 - (2*T.*alpha*gamma*rho.^(T - 1))./(rho - 1)))./(rho - 1)^2 + (24*alpha^2*gamma^2.*(T/2 - (rho.^T - 1)./(2*(rho - 1)))*(alpha + w))./(rho - 1)^3 + (12*T.*alpha^2*gamma^2*h*rho.^(T - 1))./(rho - 1)^2 - (6*T.*alpha*gamma*h.*rho.^(T - 1))./(rho - 1) - (6*T.*alpha^2*gamma^2*h*rho.^(T - 2).*(T - 1))./(rho - 1);
    %c3 = (3.*alpha.*h.*(rho.^tau - 1))/(rho - 1)^2 - (6*alpha*(tau/2 - (rho.^tau - 1)/(2*(rho - 1)))*(alpha + w))/(rho - 1)^2 - ((alpha + w)*((12*alpha^2*gamma^2*(rho.^tau - 1))/(rho - 1)^3 - (3*alpha*(rho.^tau - 1))/(rho - 1)^2 + (3*tau.*alpha.*rho.^(tau - 1))/(rho - 1) - (6*alpha*gamma*(rho.^tau - 1))./(rho - 1)^2 + (6*tau.*alpha*gamma.*rho.^(tau - 1))./(rho - 1) - (12*tau.*alpha^2*gamma^2*rho.^(tau - 1))./(rho - 1)^2 + (6*tau.*alpha^2*gamma^2.*rho.^(tau - 2).*(tau - 1))./(rho - 1)))./(rho - 1) - (12*alpha^2*gamma^2*h*(rho.^tau - 1))./(rho - 1)^3 - (3*tau.*alpha*h.*rho.^(tau - 1))./(rho - 1) + (6*alpha*gamma*h*(rho.^tau - 1))./(rho - 1)^2 - (6*alpha*gamma*(alpha + w)*(tau - (rho.^tau - 1)./(rho - 1) + (2*alpha*gamma*(rho.^tau - 1))./(rho - 1)^2 - (2*tau.*alpha*gamma*rho.^(tau - 1))./(rho - 1)))./(rho - 1)^2 + (24*alpha^2*gamma^2.*(tau/2 - (rho.^tau - 1)./(2*(rho - 1)))*(alpha + w))./(rho - 1)^3 + (12*tau.*alpha^2*gamma^2*h*rho.^(tau - 1))./(rho - 1)^2 - (6*tau.*alpha*gamma*h.*rho.^(tau - 1))./(rho - 1) - (6*tau.*alpha^2*gamma^2*h*rho.^(tau - 2).*(tau - 1))./(rho - 1);
    k3 = (3.*alpha.*h.*(rho.^tau - 1))/(rho - 1)^2 - (6.*alpha.*(tau/2 - (rho.^tau - 1)/(2.*(rho - 1))).*(alpha + w))/(rho - 1)^2 - ((alpha + w).*((12.*alpha^2.*gamma^2.*(rho.^tau - 1))/(rho - 1)^3 - (3.*alpha.*(rho.^tau - 1))/(rho - 1)^2 + (3.*tau.*alpha.*rho.^(tau - 1))/(rho - 1) - (6.*alpha.*gamma.*(rho.^tau - 1))./(rho - 1)^2 + (6.*tau.*alpha.*gamma.*rho.^(tau - 1))./(rho - 1) - (12.*tau.*alpha^2.*gamma^2.*rho.^(tau - 1))./(rho - 1)^2 + (6.*tau.*alpha^2.*gamma^2.*rho.^(tau - 2).*(tau - 1))./(rho - 1)))./(rho - 1) - (12.*alpha^2.*gamma^2.*h.*(rho.^tau - 1))./(rho - 1)^3 - (3.*tau.*alpha.*h.*rho.^(tau - 1))./(rho - 1) + (6.*alpha.*gamma.*h.*(rho.^tau - 1))./(rho - 1)^2 - (6.*alpha.*gamma.*(alpha + w).*(tau - (rho.^tau - 1)./(rho - 1) + (2.*alpha.*gamma.*(rho.^tau - 1))./(rho - 1)^2 - (2.*tau.*alpha.*gamma.*rho.^(tau - 1))./(rho - 1)))./(rho - 1)^2 + (24.*alpha^2.*gamma^2.*(tau/2 - (rho.^tau - 1)./(2.*(rho - 1))).*(alpha + w))./(rho - 1)^3 + (12.*tau.*alpha^2.*gamma^2.*h.*rho.^(tau - 1))./(rho - 1)^2 - (6.*tau.*alpha.*gamma.*h.*rho.^(tau - 1))./(rho - 1) - (6.*tau.*alpha^2.*gamma^2.*h.*rho.^(tau - 2).*(tau - 1))./(rho - 1);
    c3 = k3.*(sigma - d2)/(6*sigma^2);
    
    k4 = ((alpha + w).*((48.*alpha.^2.*gamma.*(rho.^tau - 1))./(rho - 1).^3 - (12.*alpha.*(rho.^tau - 1))./(rho - 1).^2 + (48.*alpha.^2.*gamma.^2.*(rho.^tau - 1))./(rho - 1).^3 - (96.*alpha.^3.*gamma.^3.*(rho.^tau - 1))./(rho - 1).^4 + (12.*tau.*alpha.*rho.^(tau - 1))./(rho - 1) - (48.*tau.*alpha.^2.*gamma.*rho.^(tau - 1))./(rho - 1).^2 - (48.*tau.*alpha.^2.*gamma.^2.*rho.^(tau - 1))./(rho - 1).^2 + (96.*tau.*alpha.^3.*gamma.^3.*rho.^(tau - 1))./(rho - 1).^3 + (24.*tau.*alpha.^2.*gamma.*rho.^(tau - 2).*(tau - 1))./(rho - 1) + (24.*tau.*alpha.^2.*gamma.^2.*rho.^(tau - 2).*(tau - 1))./(rho - 1) - (48.*tau.*alpha.^3.*gamma.^3.*rho.^(tau - 2).*(tau - 1))./(rho - 1).^2 + (16.*tau.*alpha.^3.*gamma.^3.*rho.^(tau - 3).*(tau - 1).*(tau - 2))./(rho - 1)))./(rho - 1) + (12.*alpha.*(alpha + w).*(tau - (rho.^tau - 1)./(rho - 1) + (2.*alpha.*gamma.*(rho.^tau - 1))./(rho - 1).^2 - (2.*tau.*alpha.*gamma.*rho.^(tau - 1))./(rho - 1)))./(rho - 1).^2 - (12.*alpha.*h.*(rho.^tau - 1))./(rho - 1).^2 + (48.*alpha.^2.*gamma.^2.*h.*(rho.^tau - 1))./(rho - 1).^3 - (96.*alpha.^3.*gamma.^3.*h.*(rho.^tau - 1))./(rho - 1).^4 - (96.*alpha.^2.*gamma.*(tau./2 - (rho.^tau - 1)./(2.*(rho - 1))).*(alpha + w))./(rho - 1).^3 - (48.*alpha.^2.*gamma.^2.*(alpha + w).*(tau - (rho.^tau - 1)./(rho - 1) + (2.*alpha.*gamma.*(rho.^tau - 1))./(rho - 1).^2 - (2.*tau.*alpha.*gamma.*rho.^(tau - 1))./(rho - 1)))./(rho - 1).^3 + (12.*tau.*alpha.*h.*rho.^(tau - 1))./(rho - 1) + (192.*alpha.^3.*gamma.^3.*(tau./2 - (rho.^tau - 1)./(2.*(rho - 1))).*(alpha + w))./(rho - 1).^4 - (8.*alpha.*gamma.*(alpha + w).*((12.*alpha.^2.*gamma.^2.*(rho.^tau - 1))./(rho - 1).^3 - (3.*alpha.*(rho.^tau - 1))./(rho - 1).^2 + (3.*tau.*alpha.*rho.^(tau - 1))./(rho - 1) - (6.*alpha.*gamma.*(rho.^tau - 1))./(rho - 1).^2 + (6.*tau.*alpha.*gamma.*rho.^(tau - 1))./(rho - 1) - (12.*tau.*alpha.^2.*gamma.^2.*rho.^(tau - 1))./(rho - 1).^2 + (6.*tau.*alpha.^2.*gamma.^2.*rho.^(tau - 2).*(tau - 1))./(rho - 1)))./(rho - 1).^2 + (48.*alpha.^2.*gamma.*h.*(rho.^tau - 1))./(rho - 1).^3 - (48.*tau.*alpha.^2.*gamma.^2.*h.*rho.^(tau - 1))./(rho - 1).^2 + (96.*tau.*alpha.^3.*gamma.^3.*h.*rho.^(tau - 1))./(rho - 1).^3 - (48.*tau.*alpha.^2.*gamma.*h.*rho.^(tau - 1))./(rho - 1).^2 + (24.*tau.*alpha.^2.*gamma.^2.*h.*rho.^(tau - 2).*(tau - 1))./(rho - 1) - (48.*tau.*alpha.^3.*gamma.^3.*h.*rho.^(tau - 2).*(tau - 1))./(rho - 1).^2 + (24.*tau.*alpha.^2.*gamma.*h.*rho.^(tau - 2).*(tau - 1))./(rho - 1) + (16.*tau.*alpha.^3.*gamma.^3.*h.*rho.^(tau - 3).*(tau - 1).*(tau - 2))./(rho - 1);
    c4 = k4/24 .* (d1.^2-1-3*sigma*d2)/sigma^3;
    
    prices = S.*normcdf(d1) - exp(-r.*tau).*K.*normcdf(d2) + S.*normpdf(d1).*exp(mu+sigma^2/2 - r*tau).*(c3+c4);
    
%     for(ind=1:nos)
%         prices(ind) =  Ralpha(F(ind), K(ind), alphas(ind))+1/pi*integral(@(x) phi(x, K(ind), alphas(ind), F(ind), kappa, theta, rho, sigma, tau(ind), v0) , 0, Inf);
%     end

    

end

