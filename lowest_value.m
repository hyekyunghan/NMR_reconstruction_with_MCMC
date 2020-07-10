clear;

phis = [0,30,60,90,300,330];
item = load(['hybrid', num2str(phis(1)) ,'.txt']);

N = size(item,1);
P = length(phis);

Nh = N/2;
xp = grid_norm(1:N,N); % make as -255.5~255.5

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
% proj = proj./mean(proj,1)*max(mean(proj,1));

%% Rotation initialize

mat = ones(N,N)*100;
mat(sqrt(xd.^2+yd.^2)>=Nh-1) = 0;

%  surf(xd(ds_idx,ds_idx), yd(ds_idx,ds_idx), mat(ds_idx,ds_idx))


%% filtering 
surf(xd(ds_idx,ds_idx), yd(ds_idx,ds_idx), mat(ds_idx,ds_idx))


for i=1:length(phis)
   
    % rotation matrix
    H = [cos(phis(i)) -sin(phis(i)); sin(phis(i)) cos(phis(i))]; 
    xr = x*H;
    
    xr_scale = zeros(length(xr),1);
    
    % scale x_coordinate after rotation
    for j=1:length(xr)
        tmp = abs(xp-xr(j,1));
        [~, min_idx] = min(tmp);
        xr_scale(j) = xp(min_idx);
    end
    xr_scale = [xr_scale, xr(:,2)];
    
    for k=1:length(proj(:,i)) % xp -255.5~255.5 / proj 1~512
        mat(xr_scale(:,1)==xp(k)) = min(mat(xr_scale(:,1)==xp(k)), proj(k,i));
    end
end
    
surf(xd(ds_idx,ds_idx), yd(ds_idx,ds_idx), mat(ds_idx,ds_idx))





