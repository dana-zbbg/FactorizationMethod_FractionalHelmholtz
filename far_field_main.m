% Plot of Far Field pattern in matrix form for Ninc incident directions
% u^infty( hat x, d)

%exponent
s = 0.6;
m = floor(1/(2*s));
%Wave number k
k = 1;
%limiting absorption parameter
eps = 1e-3;

%Number of points per row = Number of points per column
Nx = 401;
xmax = 4;
h = (2*xmax)/(Nx-1);
%Grid 
X = linspace(-xmax,xmax,Nx);
%Number of incident angles
Ninc = 36;

%Boomrang
% c = 1.6;   % small c → narrow in y, moderately sized in x
% a = 0.8; 
% inhom = @(x,y) 10*(((c*x + 2*a/c*y.^2).^2+y.^2 ) <= c^2 & (x+2.5).^2 + y.^2 <=4);

% bunny 2.0
% inhom = @(x,y) 10*((x.^2*16/9 + (y+0.5).^2*4 <= 1) |...
%     ( ((x - 0.6)*0.5 + (y-0.3)*0.5*sqrt(3)).^2*2 ...
%     +(-(x - 0.6)*0.5*sqrt(3) + (y-0.3)*0.5).^2*16 <= 1) | ...
%     ( (-(x + 0.6)*0.5 + (y-0.3)*0.5*sqrt(3)).^2*2 + ...
%     (-(x + 0.6)*0.5*sqrt(3) - (y-0.3)*0.5).^2*16 <= 1));

% cercle+ rectangle
% inhom = @(x,y) 20 * ( (x-1).^2 + (y-1).^2 <= 0.25) + ...
%    10* ((-1 <= x)& (x <= -0.5) & (-1 <= y)& (y<=0));


% kite
c= 1.6;
a= 0.9;
inhom = @(x,y) 10*(((c*x + 2*a/c*y.^2).^2+y.^2 ) <= c^2);


%% inhomogeneity
F = inhom(X, X');
mask = abs(F) > 1e-16;
[y_support, x_support] = find(mask);
nb_support = numel(x_support);

inc = @(x1,y1,x2,y2) exp(1i*k*(x1*x2+y1*y2));
%% Lippmann Schwinger Matrix
tic
[F,M,Gsh] = LippmannSchwingerMatrix_assemble(s,k,h,Nx,xmax,X,inhom, Ninc, nb_support, x_support,y_support);
t0 = toc;
disp(['Elapsed time for LippmannSchwinger : ', num2str(t0), ' seconds ']);
tic
luM = decomposition(M, "lu");
%Minv= inv(M);
time_M_LU = toc;
disp(['Elapsed time for LU decomposition of M = LippmannSchwinger Matrix : ', num2str(time_M_LU), ' seconds ']);

%% Angles of the incident field
angle_polar = linspace(0,2*pi*(1-1/Ninc), Ninc);
angles_cos = cos(angle_polar);
angles_sin = sin(angle_polar);

%% Far Fields Initialization
FarFieldFF = zeros(Ninc, Ninc);

%% for each incident direction compute Far Field
for theta= 1:Ninc
     % initializing RHS and solution 
    tic
    u_inc = zeros(nb_support,1);
    for p1=1:nb_support
        u_inc(p1) = inc(X(x_support(p1)), X(y_support(p1)), angles_cos(theta) , angles_sin(theta));
    end
    u_tot = luM\u_inc;
    FarFieldFF(:, theta) = F*u_tot;
end

save('FF_data.mat', 'k', 's','inhom', "FarFieldFF")




hFig = figure('Visible','off');
pcolor(real(FarFieldFF), EdgeColor="none");
colorbar;
xlabel("angle")
ylabel("Far Field")
title("k=20, far field for a small inhomogeneity, s = "+ num2str(s));
filename = sprintf('farfield_matrix.png');
exportgraphics(hFig, filename, 'Resolution', 300);
clf(hFig);
close(hFig);