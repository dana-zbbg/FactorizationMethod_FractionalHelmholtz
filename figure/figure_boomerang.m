
% Plot reconstruction of the 'boomerang' inhomogeneity, for
% different values of k and Ninc

%% True scatterer
c = 1.6;   
a = 0.8; 
inhom = @(x,y) 10*(((c*x + 2*a/c*y.^2).^2+y.^2 ) <= c^2 & (x+2.5).^2 + y.^2 <=4);

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
p1 = load('FM_boomerang_s=0,6_Ninc=36_k=2.mat');
norm1 = max(abs(p1.W(:)));
contourf(x, x, p1.W./norm1, 30, 'LineColor', 'none');
colormap(ax2,1-gray)
title('Factorization Method k=2, Ninc=36');


ax3 = subplot(2, 2, 3);  
p2 = load('FM_boomerang_s=0,6_Ninc=36_k=5.mat');
norm2 = max(abs(p2.W(:)));
contourf(x, x, p2.W./norm2, 30, 'LineColor', 'none');
colormap(ax3,1-gray)
title('Factorization Method k=5, Ninc=36');


ax4 = subplot(2, 2, 4); 
p5 = load('FM_boomerang_s=0,6_Ninc=72_k=5.mat');
norm5 = max(abs(p5.W(:)));
contourf(x, x, p5.W./norm5, 30, 'LineColor', 'none');
colormap(ax4,1-gray)
title('Factorization Method k=5, Ninc=72')

cb = colorbar(ax4);
cb.Position(1) = ax4.Position(1) + ax4.Position(3) + 0.01;
cb.Position(2) = ax3.Position(2);
cb.Position(4) = ax2.Position(2) + ax2.Position(4) - ax3.Position(2);


% filename = sprintf('fig_boomerang.png');
% exportgraphics(hFig, filename, 'Resolution', 300);