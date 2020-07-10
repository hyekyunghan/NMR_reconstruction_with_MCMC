function [y_hat] = func_yhat(xp,mu,A,phis)
%FUNC_YHAT Summary of this function goes here
%   Detailed explanation goes here
    H = zeros(2,length(phis));
    N = length(xp);
    K = length(mu);

    for i=1:length(phis)
        H(:,i) = [cos(phis(i)) ; sin(phis(i))];
    end
    mu = mu*H;

    y_hat = zeros(length(phis),N); 

    for j=1:length(phis) %5
        for i=1:K
            y_hat(j,:) = y_hat(j,:) + A(i)*normpdf(xp,mu(i,j),2);
        end
    end
end

