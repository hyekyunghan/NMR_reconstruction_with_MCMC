function lkhs = likelihd2(xp,y,mu,A,sigma_e)
%LIKELIHD Summary of this function goes here
%   Detailed explanation goes here

phis = [0,30,60,90,300,330];

% y_hat = zeros(length(phis),N); % vector
% lkh = 0;
% lkhs = [];
% H = zeros(2,length(phis));


%% creating rotation matrix
% for i=1:length(phis)
%     H(:,i) = [cos(phis(i)) ; sin(phis(i))];
% end

%%
y_hat = func_yhat(xp,mu,A,phis);
% mu = mu*H;

% y_hat = zeros(length(phis),N); 
% 
% 
% 
% for j=1:length(phis) %5
%     for i=1:K
%         y_hat(j,:) = y_hat(j,:) + A(i)*normpdf(xp,mu(i,j),2);
%     end
%     
% end

dif = normpdf(y-y_hat,0,sigma_e);
lkhs = sum(log(dif(:)));

% keyboard
% lkh = normpdf(y(j)-y_hat,0,sigma_e);
% lkh = sum(log(lkh));
% lkhs(j) = lkh;
% lkhs = sum(lkhs);