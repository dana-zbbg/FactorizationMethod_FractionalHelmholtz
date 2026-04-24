
% Plot reconstruction of the 'cercle+rectangle' inhomogeneity, for
% different values of s

%% True scatterer : cercle+rectangle
inhom = @(x,y) 20 * ( (x-1).^2 + (y-1).^2 <= 0.25) + ...
   10* ((-1 <= x)& (x <= -0.5) & (-1 <= y)& (y<=0));
%% Parameters
Nx = 401;
xmax = 5;
x = linspace(-xmax,xmax, Nx);
M = inhom(x,x');
Mn = M/max(M(:));

%% Figure
hFig = figure;

ax1 = subplot(2, 2, 1);  
imagesc(x,x,Mn);
colormap(gca, 1-gray);   
title('True Scatterer');

ax2 = subplot(2, 2, 2);  
p1 = load('FM_cercle+rectangle_s=0,2_Ninc=72_k=5.mat');
norm1 = max(abs(p1.W(:)));
contourf(x, x, flipud(p1.W./norm1), 30, 'LineColor', 'none');
colormap(ax2,1-gray)
title('Factorization Method s=0,2');


ax3 = subplot(2, 2, 3);  
p2 = load('FM_cercle+rectangle_s=0,7_Ninc=72_k=5.mat');
norm2 = max(abs(p2.W(:)));
contourf(x, x, flipud(p2.W./norm2), 30, 'LineColor', 'none');
colormap(ax3,1-gray)
title('Factorization Method s=0,7');


ax4 = subplot(2, 2, 4); % third row right
p3 = load('FM_cercle+rectangle_s=1_Ninc=72_k=5.mat');
norm5 = max(abs(p3.W(:)));
contourf(x, x, flipud(p3.W./norm5), 30, 'LineColor', 'none');
colormap(ax4,1-gray)
title('Factorization Method s=1 (Helmholtz)')

cb = colorbar(ax4);
cb.Position(1) = ax4.Position(1) + ax4.Position(3) + 0.01;
cb.Position(2) = ax3.Position(2);
cb.Position(4) = ax2.Position(2) + ax2.Position(4) - ax3.Position(2);

%filename = sprintf('fig_cercle+rectangle.png');
%exportgraphics(hFig, filename, 'Resolution', 300);