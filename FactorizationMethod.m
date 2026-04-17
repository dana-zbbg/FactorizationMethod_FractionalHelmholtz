%% =========================================================
%  FACTORIZATION METHOD FOR INVERSE SCATTERING (2D)
%  ---------------------------------------------------------
%  Reconstructs the support of a scatterer from far-field
%  data using the Factorization Method (Kirsch 1998).
%
%  Problem: time-harmonic acoustic scattering
%    ((-Delta)^s - k^2s) u^s = k^2s (n-1) u
%
%  The far-field operator F has the factorization:
%       F = G S G*
%  The indicator is:  W(z) = 1 / sum_j ( |<phi_z, v_j>|^2 / |lambda_j| )
%  z is INSIDE D  <=>  phi_z is in the range of |F|^{1/2}
% =========================================================

% to use pre-computed far field data use
addpath '...\far_field_data\'


clear; close all; clc;

load("FF_data.mat")

Ninc = length(FarFieldFF.'); 
alpha   = 1e-3;       % Tikhonov regularisation (for SVD truncation)

% Sampling grid
Ngrid   = 401;        % grid points per side
L       = 5;          % half-width of sampling domain
x1      = linspace(-L, L, Ngrid);
[X, Y]  = meshgrid(x1, x1);
ZZ      = [X(:), Y(:)];   % all sampling points  (Ngrid^2 x 2)

theta   = (0:Ninc-1)' * 2*pi / Ninc;   % incident / observation angles
d       = [cos(theta), sin(theta)];   % direction vectors  (N x 2)

fprintf('Computing SVD ...\n');
[U, S_mat, V] = svd(FarFieldFF); 
sigma = diag(S_mat);          % singular values  (N x 1)

% Picard plot 
figure('Name','Picard plot');
semilogy(sigma,'b.-','LineWidth',1.5,'MarkerSize',10);
xlabel('Index j'); ylabel('|\sigma_j|');
title('Singular values of F'); grid on;

fprintf('Computing indicator on %d x %d grid ...\n', Ngrid, Ngrid);
% Indicator accumulator
sigma_abs  = abs(sigma);
sigma_tik  = sigma_abs ./ (sigma_abs.^2 + alpha^2);   % Tikhonov weights
%M = 25;   % truncate at the cliff
%sigma_tik = zeros(size(sigma_abs));
%sigma_tik(1:M) = 1 ./ sigma_abs(1:M);   % zero out the rest

W = zeros(Ngrid, Ngrid);
for ix=1:Ngrid
    ax=1i*k*x1(ix)*cos(theta);
    for iy=1:Ngrid
        sc=V'*exp(-(ax+1i*k*x1(iy)*sin(theta))); % V'
        GUP=conj(sc).*sc.*sigma_tik;
        W(iy,ix)=1/sqrt(sum(GUP));
    end
end   

Fig1 = figure('Name','Factorization Method','Position',[100 100 1200 450]);
ncount =5;
t = tiledlayout(2,3);
%--- Left: Indicator function
%% ----TILE 1 
ax1 = nexttile;
contourf(x1, x1, W,ncount);
title(ax1,'Result of Factorization Method')
colormap(ax1,1-gray)
colorbar(ax1)
axis(ax1,'square')
%% --- TILE 2 ---
ax2 = nexttile;
[X1,X2] = meshgrid(x1,x1);
M = inhom(X1,X2);
Mn = M/max(M(:));
set(ax2,'YDir','normal')
colormap(ax2,[ones(256,1), linspace(1,0,256)', linspace(1,0,256)'])
caxis(ax2,[0 1])
title('True Scatterer')
colorbar(ax2)

hold(ax2,'on')
h = imagesc(ax2,x1,x1,ones(size(Mn)));
set(h,'AlphaData',0.6*Mn)
axis(ax2,'square')
%% ----TILE 3
ax3 = nexttile;
contourf(ax3,x1,x1,W,ncount)
title(ax3,'Factorization Method Result and Scatterer')
colormap(ax3,1-gray)
colorbar(ax3)
axis(ax3,'square')

hold(ax3,'on')
RGB = cat(3,ones(size(M)),zeros(size(M)),zeros(size(M)));
h = imagesc(ax3,x1,x1,RGB);
set(ax3,'YDir','normal')
set(h,'AlphaData',0.6*Mn)

sgtitle(sprintf('Factorization Method  (k = %g,  Ninc = %d, s=%g)', ...
    k, Ninc,s), 'FontSize', 14);

Fig2 = figure;
contourf(x1, x1, W, 30, 'LineColor', 'none');
colormap hot; colorbar;
set(gca, 'YDir', 'normal');  % or axis xy
sgtitle(sprintf('Factorization Method  (s = %g)', ...
    s), 'FontSize', 14);

save('FM.mat', 'W')

% filename1 = sprintf('fig_kite.png');
% exportgraphics(Fig1, filename1, 'Resolution', 300);
% filename2 = sprintf('fig_kite_hot.png');
% exportgraphics(Fig2, filename2, 'Resolution', 300);


