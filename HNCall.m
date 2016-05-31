function [ call, i1, i2 ] = HNCall( S, X, T, sigma, rf, alpha, beta, gamma, lambda, omega)
%HNCALL Summary of this function goes here
%   Detailed explanation goes here
    lambda_rn = -.5;
    %lambda_rn = lambda;
    gamma_rn = gamma + lambda + .5;
    %gamma_rn = gamma;

    i1 = (exp(-rf*T)/pi) * quad (@Integral1, eps, 100) ;
    i2 = (0.5 + quad (@Integral2, eps, 100)/pi) ;
    
    call = .5*S + (exp(-rf*T)/pi) * quad(@Integral1, eps, 100) - X*exp(-rf*T)*(0.5 + quad(@Integral2, eps, 100)/pi);

    function f1 = Integral1(phi)
        f1 = real((X.^(-1i*phi) .* genFunc(1i*phi+1)) ./ (1i * phi));
    end

    function f2 = Integral2(phi)
        f2 = real((X.^(-1i*phi) .* genFunc(1i*phi)) ./ (1i * phi));
    end

    function f = genFunc(phi)
        phi = phi';
        
        A(:,T) = phi * 0;
        B(:,T) = phi * 0;
        for i=T-1:-1:1
            A(:,i) = A(:,i+1) + phi.*rf + B(:,i+1).*omega - .5 * log(1- 2 * alpha * B(:,i+1));
            B(:,i) = phi.*(lambda_rn + gamma_rn) - .5 * gamma_rn^2 + beta .* B(:, i+1) + .5*(phi-gamma_rn).^2./(1-2*alpha*B(:,i+1));
        end
        
        A_ = A(:,1) + phi.*rf + B(:,1).*omega - .5 * log(1- 2 * alpha * B(:,1));
        B_ = phi.*(lambda_rn + gamma_rn) - .5 * gamma_rn^2 + beta .* B(:, 1) + .5*(phi-gamma_rn).^2./(1-2.*alpha.*B(:,1));
        
        f = S.^phi .* exp(A_+sigma.*B_);
        f = f';
    end

end
