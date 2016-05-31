function [ call ] = HNCall( S_, X_, T_, sigma_, rf_, model)
%HNCALL Summary of this function goes here
%   Detailed explanation goes here
    gamma = model.gamma;
    lambda = model.lambda;
    alpha = model.alpha;
    beta = model.beta;
    omega = model.omega;
    lambda_rn = -.5;
    %lambda_rn = lambda;
    gamma_rn = gamma + lambda + .5;
    %gamma_rn = gamma;

    call = .5*S_ + (exp(-rf_*T_)/pi) * quad(@Integral1, eps, 100) - X_*exp(-rf_*T_)*(0.5 + quad(@Integral2, eps, 100)/pi);

    function f1 = Integral1(phi)
        f1 = real((X_.^(-1i*phi) .* genFunc(1i*phi+1)) ./ (1i * phi));
    end

    function f2 = Integral2(phi)
        f2 = real((X_.^(-1i*phi) .* genFunc(1i*phi)) ./ (1i * phi));
    end

    function f = genFunc(phi)
        phi = phi';
        
        A(:,T_) = phi * 0;
        B(:,T_) = phi * 0;
        for i=T_-1:-1:1
            A(:,i) = A(:,i+1) + phi.*rf_ + B(:,i+1).*omega - .5 * log(1- 2 * alpha * B(:,i+1));
            B(:,i) = phi.*(lambda_rn + gamma_rn) - .5 * gamma_rn^2 + beta .* B(:, i+1) + .5*(phi-gamma_rn).^2./(1-2*alpha*B(:,i+1));
        end
        
        A_ = A(:,1) + phi.*rf_ + B(:,1).*omega - .5 * log(1- 2 * alpha * B(:,1));
        B_ = phi.*(lambda_rn + gamma_rn) - .5 * gamma_rn^2 + beta .* B(:, 1) + .5*(phi-gamma_rn).^2./(1-2.*alpha.*B(:,1));
        
        f = S_.^phi .* exp(A_+sigma_.*B_);
        f = f';
    end

end
