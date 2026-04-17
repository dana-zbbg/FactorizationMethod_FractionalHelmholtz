%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Lippmann Schwinger Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computation of the solution
function  [F,M,Gsh] = LippmannSchwingerMatrix_assemble(s,k,h,Nx,xmax, X, inhom, Ninc, nb_support, x_support, y_support)% xq_beg, xq_end, yq_beg, yq_end)
    %% computation of green function
    m = floor(1/(2*s));
    %Helmholtz part i/4 * H^(1)_0(k |x|)
    Ghelm_fun = @(d) s^(-1)*k^(2-2*s)*1i*0.25*besselh(0, k*d);
    if s>0.5
        F = @(R) ((R^2 - k^2 - s^(-1)*k^(2-2*s)*(R^(2*s)-k^(2*s)))*R)/((R^2-k^2)*(R^(2*s)-k^(2*s)));
    elseif s< 0.5
        if 1/(2*s) == m%floor(1/(2*s))
            F = @(R) (k*(R^2 - k^2) + k^(2-2*s)*(R^(2*s)-k^(2*s))*(R-k-s^(-1)*R)) ...
                /((R^2-k^2)*(R^(2*s)-k^(2*s)));
        else 
            F = @(R) (R^(1-2*s*m)*(k^(2*s*m)*(R^2 - k^2) - s^(-1)*k^(2-2*s)*R^(2*s*m) ...
                *(R^(2*s)-k^(2*s))))/((R^2-k^2)*(R^(2*s)-k^(2*s)));
        end
    end
    dist = linspace(0,2*sqrt(2)*xmax, 2*Nx);
    h0 = 2*sqrt(2)*xmax/(2*Nx-1);
    tic
    max_period = 2*pi/(2*sqrt(2)*xmax);
    Nb_points_per_period = 10;
    rmax = max(1e2, floor(h^(-2/(2*s*(m+1)-0.5))));
    Nb_periods = floor(rmax/max_period);
    Nr = Nb_periods*Nb_points_per_period;
    disp(['Nr = ', num2str(Nr), ' rmax = ', num2str(rmax)]);
    r = linspace(0,rmax, Nr);%val = linspace(0,1e3,1e6);

    
    Ghelm = arrayfun(Ghelm_fun, dist);
    if s==0.5
        Gsh = Ghelm;
    elseif s==1
        Gsh = Ghelm;
    else
        FX = arrayfun(F,r);
        Gs = hankel_transform(FX,r,dist);
        Gsh = Gs + Ghelm;
    end

    if s<=0.5
        for j=0:(m-1)
            power = @(x) k^(2*s*j)*gamma(1-s*(j+1))/(4^(s*(j+1))*...
            pi*gamma(s*(j+1))*x^(2-2*s*(j+1)));
            Gsh = Gsh + arrayfun(power, dist);
        end
        if 1/(2*s) == m
            struveK = @(x) -k^(2-2*s)*0.25*StruveH0Y0(k*x);
            Gsh = Gsh + arrayfun(struveK, dist);
        end
    end

    %Integral on ball i,j of the Kernel
    %near 0, i/4 * H^(1)_0(k |x|) \sim i/4 * 2i/\pi *ln(k|x|) etc
    IGhelm = (1./16)*s^(-1)*k^(2-2*s)*(1 - 2*log(k*h./2)) -(1./8)*s^(-1)*k^(2-2*s)*0.5772;
    %IGhelm = -(1./16)*s^(-1)*k^(2-2*s)*(-1 + 2*log(k*h./2));
    IGStruve = -(1./16)*k^(2-2*s)*(2*log(h*k/2) - 1)  - k^(2-2*s)*0.57721/2;
    IGintegral = h^(2*s*(m+1)-2)*k^(2*s*m)*4^(-2*s*(m+1))*gamma(1-s*(m+1))/gamma(1+s*(m+1));%...
           % + h^(4*s-2)*k^(2*s)*2^(-4*s)*gamma(1-2*s)/gamma(1+2*s);%test for more accuracy;
    %(singularity from kernels done after)
    % s<0.5 integer case : struve and integral contribution 
    if 1/(2*s)==m
        if s==0.5
            IGs = IGStruve;
        else
            IGs = IGintegral + IGStruve;
        end
    % otherwise 
    elseif s==1
        IGs = 0;
    else
        IGs = IGintegral;
    end
    
    Gsh(1) = IGhelm+IGs;
    disp(num2str(Gsh(1)));
    if s <= 0.5
        for j=0:(m-1)
            Gsh(1) = Gsh(1) + k^(2*s*j)*gamma(1-s*(j+1))/(4^(s*(j+1))*...
            s*(j+1)*gamma(s*(j+1)))*(0.5*h)^(2*s*(j+1));
        end
    end
    time_green = toc;
    disp(['Elapsed time for green function : ', num2str(time_green), ' seconds ']);

%% computation of scattered field
    tic
    %[nb_support, x_support, y_support] = support_xy(inhom, X);
    %[support] = obtain_support(nb_support, x_support, y_support,Nx);
    %number of coefficients in the sparse matrix
    disp([" Size Matrix : ", nb_support^2]);
    if nb_support^2 > 8*1e9
        disp("Error : size matrix too large");
        return
    end
    M = zeros(nb_support, nb_support);
    tic
    [xi_ind, yi_ind, xj_ind, yj_ind, xi,yi,xj, yj]=deal(1);
    indexG = 0;
    for p1=1:nb_support
        for p2=1:nb_support
            xi_ind = x_support(p1);%support(p1,2);
            yi_ind = y_support(p1);%support(p1,3);
            xj_ind = x_support(p2);%mod(support(p2)-1,Nx)+1;
            yj_ind = y_support(p2);%floor((support(p2) - xj_ind)/Nx) + 1;
            
            xi = X(xi_ind);
            yi = X(yi_ind);
            xj = X(xj_ind);
            yj = X(yj_ind);
            %dist between (xi,yi) and (xj, yj)
            r = sqrt((xi-xj)^2+(yi-yj)^2); %at most sqrt(2)*(2Nx-1)
            indexG = floor(r/h0)+1;
            M(p1,p2) = -k^(2*s)*Gsh(indexG)*h^2*inhom(xj,yj)*h^2;
            if p1 == p2
               M(p1,p2) = M(p1,p2) + 1;
            end
        end
    end
    time_M = toc; 
    disp(['Elapsed time for I-M : ', num2str(time_M), ' seconds ']);
    F = zeros(Ninc, nb_support);
    inc = @(a1,a2,y1,y2) exp(-1i*k*(a1*y1+a2*y2));
    angles = linspace(0,2*pi*(1-1/Ninc), Ninc);
    for a=1:Ninc
        a1 = cos(angles(a));
        a2 = sin(angles(a));
        for p2=1:nb_support
            xj_ind = x_support(p2);%support(p2,2);%mod(support(p2)-1,Nx)+1;
            yj_ind = y_support(p2);%support(p2,3);%floor((support(p2) - xj_ind)/Nx) + 1;
            
            xj = X(xj_ind);
            yj = X(yj_ind);
            F(a,p2) = k^2/s*inc(a1,a2,xj,yj)*h^2*inhom(xj,yj);
        end
    end

end
