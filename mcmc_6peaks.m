clear;

phis = [0,30,60,90,300,330];
item = load(['hybrid', num2str(phis(1)) ,'.txt']);

N = size(item,1);
P = length(phis);

Nh = N/2;
% Nd = (1:N)-(N+1)/2;
xp = grid_norm(1:N,N);

ds_idx = 1:2:N; % Downsampling idxxd

[xd,yd] = meshgrid(xp,xp);
proj = zeros(N,length(phis));

for i=1:length(phis)
    phi = phis(i);
    item = load(['hybrid', num2str(phi) ,'.txt']);
    proj(:,i) = item(:,2);
end

%%
x = [xd(:), yd(:)];
phis = deg2rad(phis);
% [V,I] = max(mean(proj,1));
proj = proj./mean(proj,1)*max(mean(proj,1));

%% Rotation initialize

mat = ones(N,N)*100;
mat(sqrt(xd.^2+yd.^2)>=Nh-1) = 0;

%  surf(xd(ds_idx,ds_idx), yd(ds_idx,ds_idx), mat(ds_idx,ds_idx))

%% set prior
% mu ~ uniform
% A ~ normal
% sigma ~ Inverse Gamma
% mu_prior = unifpdf(xp,min(xp),max(xp));
% A_prior = normpdf(xp,150,2);
% sigma_e_prior = 1./gampdf(xp,10,2);

%% initialising parameters

K = 5;
mu = unifrnd(min(xp),max(xp),[K,2]);
A = unifrnd(0,300,[1,K]);
sigma_e = 20;

%% creating rotation matrix
H=[];
for i=1:length(phis)
    h = [cos(phis(i)) ; sin(phis(i))];
    H = [H h];
end


%% my model 
y = proj';


%% posterior 

rng(3);
mu_sig = 20;
A_sig = 15;
sigma_e_sig = 5;
iteration = 1000;
hist=zeros(16,length(iteration));
hist_pstr=zeros(1,length(iteration));
for i=1:iteration
    for k=1:5
        
        % mu update
%         for h=1:length(y(:,1))
        mu_prior = unifpdf(mu(k,:),min(xp),max(xp));
        mu_star = mu;

        mu_star(k,:) = normrnd(mu(k,:),mu_sig,[1,2]);
        mu_star_prior = unifpdf(mu_star(K),min(xp),max(xp));
        
        
        w = likelihd2(xp,y,mu_star,A,sigma_e)+log(mu_star_prior)...
            -likelihd2(xp,y,mu,A,sigma_e)-log(mu_prior);
        
        a = min(log(1),w);
        u = rand;

        if log(u) < a
            mu = mu_star;
        end
        
        % update A
        A_star = A;
        A_star(k) = normrnd(A(k),A_sig);
        
        A_prior = normpdf(A_star(k),A(k),A_sig);
        A_star_prior = normpdf(A(k),A_star(k),A_sig);

        w = likelihd2(xp,y,mu,A_star,sigma_e)+log(A_star_prior)...
            -likelihd2(xp,y,mu,A,sigma_e)-log(A_prior);
        a = min(log(1),w);
        u = rand;

        if log(u) < a
            A = A_star;
        end
    end
    
    % update sigma_e
    sigma_e_star = normrnd(sigma_e,sigma_e_sig);
    sigma_e_prior = 1./gampdf(sigma_e,10,2);
    sigma_e_star_prior = 1./gampdf(sigma_e_star,10,2);

    w = likelihd2(xp,y,mu,A,sigma_e_star)+log(sigma_e_star_prior)...
        -likelihd2(xp,y,mu,A,sigma_e)-log(sigma_e_prior);
    a = min(log(1),w);
    u = rand;

    if log(u) < a
        sigma_e = sigma_e_star;
    end

    y_hat = func_yhat(xp,mu,A,phis);
 
    
    hist(:,i) = [mu(:,1)' mu(:,2)' A sigma_e]';
    
    mu_hist=0;
    A_hist=0;
    for t=1:5
        mu_prior_hist = unifpdf(mu(t,:),min(xp),max(xp));
        mu_prior_hist = log(mu_prior_hist);
        mu_hist = mu_hist + mu_prior_hist;
        
        A_prior_hist = normpdf(A_star(t),A(t),A_sig);
        A_prior_hist = log(A_prior_hist);
        A_hist = A_hist + A_prior_hist;
    end
     
    hist_pstr(i) = likelihd2(xp,y,mu,A,sigma_e) + mu_hist(1) + ...
        A_hist + log(1./gampdf(sigma_e,10,2));

    figure(1); clf

    for i_plot = 1:length(phis)
        subplot(3,2,i_plot);
        plot(y(i_plot,:)); hold on
        plot(y_hat(i_plot,:))
    end

    pause(0.001);
    
    
    
end

% make history
% to find the optimal parameter
% where the posterior is the highest value
[M, I] = max(hist_pstr);
theta = hist(:,I);


mu_opt = zeros(5,2)
A_opt = zeros(5,1)
for i=1:5
    mu_opt(i,:) = [hist(i,I) hist(i+5,I)];
end
for i=1:5
    A_opt(i) = hist(i+10,I);
end
sigma_e_opt = hist(16,I);

%% create image looking from the above
% White dots are the peaks we've been looking for
mat2 = zeros(length(xp)*length(xp),1);
for k = 1:K
    mat2 = mat2+ mvnpdf(x,mu(k,:),[sigma_e 0;0 sigma_e])*A(k);
end

mat2 = reshape(mat2, size(mat));
figure(2);
imshow(mat2)


