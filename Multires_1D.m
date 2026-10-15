%Multi-resolution WENO implementation in 1D (constrained) 
    %%% Run nonlinear = 0 for order or nonlinear = 1 for behavior around shock

clear all;close all;format short e;

Method = "multi-resolution WENO Scheme";
order = 3; %Set desired error order (3 or 5)
nonlinear = 0;


%r = 2; k = 3; m = 5; %
r = 1; k = 3; m = 3; %

weights = [1/11, 10/11]; %


%Uncomment the desired time
%t_final = 0.3; %Before shock
t_final = 2; %After shock

%Comment or uncomment if you want the shock for T=2.
avoid_shock = true; %(without shock)
%avoid_shock = false; %(with shock)

%Classic, q_k smoothness (1) or multi-res smoothness (2) 
smoothness = 1;

%Defining f and initial condition for Burger's 

alpha = 1/3; beta = 2/3; 
u_initial = @(x) alpha + beta*sin(x); 
xL = 0; xR = 2*pi;

%Generate points for exact solution (only do once)
x0 = linspace(pi/400,pi-pi/400,3657); 

%Range of points on x for convergence study
n_array=[20,40,80,160]; 
% n_array=[600,700,800]; 
%n_array=[50]; 

%Preallocate arrays for memory
dx = zeros(length(n_array),1); 
L2Error = zeros(length(n_array),1); L1Error = zeros(length(n_array),1); 
LinfError= zeros(length(n_array),1); err = zeros(length(n_array),1); 

%Cycle through all ranges of dx
for kk=1:length(n_array)

    %Build spatial array
    N = n_array(kk);
    dx(kk)=(xR-xL)/N;
    h=dx(kk); 
    x = (xL:h:xR)'; %We have nSteps + 1 points here
    x_avg = (x(1:end-1)+x(2:end))/2;
    
    %Preallocate arrays
    u_avg=zeros(N,1); u_ref = zeros(size(u_avg));
    
    %Initialize cell average and fill in initial conditions
    for j=1:N
        [vals,weights] = lgwt(2,x(j),x(j+1));
        %u_avg(j) = integral(u_initial, x(j), x(j+1), 'ArrayValued', true)/h;
        u_avg(j) = one_dim_gauss_quad(u_initial,vals,weights)/h;
    
    end 

    %Obtain weight matrix for coefficents
    W = zeros(k, m+1, N);
    for i = 1:N
        W(:,:,i) = condensed_constrained_weight_matrix(k,m,r,i,x);
    end
    A = zeros(m+1, m+1, N);
    for i = 1:N
        A(:,:,i) = stiffness_weight_matrix(k,m,r,i,x,W(:,:,i));
    end

    t = 0;
    while t < t_final     
        %Choose appropriate time step
        a = max(abs(u_avg));
        if order == 3
            dt = (h)/(2*a); 
        elseif order == 5
            dt = h^(5/3); %Try to get fifth order acuracy by taking smaller dt
        end

        if t_final - t < dt
             dt = t_final - t;
        end

        l= dt/(h);

        %RK3
        u_hold = L(u_avg,weights,a,l,W,A,x,smoothness);
        u_hold = 3/4*u_avg + 1/4*L(u_hold,weights,a,l,W,A,x,smoothness);
        u_next = 1/3*u_avg + 2/3*L(u_hold, weights,a,l,W,A,x,smoothness);

        u_avg = u_next;
        t = t + dt;
    end
    
    
    exact_sol = @(x) solveBurgers2(alpha, beta, t_final, x);
    for j=1:N 
       u_ref(j) = 1/h*integral(exact_sol, x(j), x(j+1), 'ArrayValued', true);
    end

    
    %Error and convergence analysis

    if avoid_shock == false
    uDif=u_ref(:)-u_avg(:);
        L1Error(kk)= (1/N)*sum(abs(uDif));
        L2Error(kk)= ((1/N)*sum(abs(uDif).^2))^(1/2);
        LinfError(kk) = max(abs(uDif));
    end
    
    if avoid_shock == true
        % Error from points far from shock
        xShock = pi + alpha * t_final;
        idx = find(abs(x_avg-xShock)>1);
    
        uDif = u_ref(idx) - u_avg(idx);
        L1Error(kk)= (1/numel(idx))*sum(abs(uDif));
        L2Error(kk)= (1/numel(idx))*sum(abs(uDif).^2)^(1/2);
        LinfError(kk) = max(abs(uDif));
    end




    %Plotting at desired time
    variable_name = "Burger's";
    if N == 160
        figure(kk)
        plot(x_avg,u_avg,'.-')
        hold on
        plot(x_avg, u_ref, 'ro')
        title(sprintf("Plot of %s from %s with number of points ", [variable_name Method]) + N)
        legend('Approximate','Exact')

        %saveas(kk, sprintf('%s/1D/constrained_multres_tFinal_%d_smoothness_%d_plot.jpg',pth,t_final,smoothness),'jpeg')
    end
end


p = order;

fig_num = 9; 

 order_plot("Burgers" ,L2Error,fig_num,Method,p,dx,"L2 Error", k)

 if t_final == 0.3
 %saveas(fig_num, sprintf('%s/1D/constrained_multres_tFinal_%d_smoothness_%d_convergence_plot.jpg',pth,t_final,smoothness),'jpeg')
 end


 if t_final == 2
 %saveas(fig_num, sprintf('%s/1D/constrained_multres_tFinal_%d_smoothness_%d_avoid_shock_%d_convergence_plot.jpg',pth,t_final,smoothness,avoid_shock),'jpeg')
 end

 order_table(dx,n_array, L2Error, "L2 Error");
 order_table(dx,n_array, L1Error, "L1 Error");
 order_table(dx,n_array, LinfError, "Linf Error");




%Functions used (for now)

function u_hold = L(u,weights,a,l,W,A,x,smoothness)
    [uphm,umhp] = wenoReconstruction(u,weights,W,A,x,smoothness);
    uphp = circshift(umhp, -1); %
    fhatp = lf(a,uphm, uphp); %
    fhatm = circshift(fhatp,1);
    u_hold = (u - l.*(fhatp - fhatm));
end

function [uphm,umhp]  = wenoReconstruction(u_avg,weights,W,A,x,smoothness)
%wI1 = weights(1); wI2 = weights(2); wI3 = weights(3);
epsilon=1e-5;
%gam_12 = weights(1); gam_22 = weights(2);
gam_12 = 1/11; gam_22 = 10/11;

N = length(u_avg);

dx = circshift(x, -1) - x;
up2 = circshift(u_avg, -2);
up1 = circshift(u_avg, -1); 
um1 = circshift(u_avg,  1);
um2 = circshift(u_avg,  2);

uphm = zeros(N,1); umhp = zeros(N,1);

for i = 1:N
    h = dx(i);

    uip2 = up2(i); uip1 = up1(i); uim1  = um1(i); uim2  = um2(i); ui = u_avg(i);
    xph = x(i + 1); xmh = x(i);
    
  
    %d = [uim2; uim1; ui; uip1; uip2 ; -ui];
    d = [uim1; ui; uip1 ; -ui];

    C = W(:,:,i)*d;
    c1 = C(1); c2 = C(2); c3 = C(3);

    q1 = @(x) ui;
    q2 = @(x) c1 + c2*x + c3*x^2;

    p1 = @(x) q1(x);
    p2 = @(x) (11/10)*q2(x) - (1/10)*p1(x);

    S_0 = (ui - uim1)^2; S_1 = (uip1 - ui)^2;

    if S_0 >= S_1
        gam_bar_01 = 1;
    else 
        gam_bar_01 = 10;
    end

    gam_bar_11 = 11 - gam_bar_01;

    gam_01 = gam_bar_01/(gam_bar_01 + gam_bar_11);
    gam_11 = 1 - gam_01;

    sigma_0 = gam_01 *(1 + abs(S_0 - S_1)/(S_0 + epsilon));
    sigma_1 = gam_11*(1 + abs(S_0 - S_1)/(S_1 + epsilon));
    sigma =   sigma_0 + sigma_1;

    beta_1 = (1/sigma^2)*(sigma_0*(ui - uim1) + sigma_1*(uip1 - ui))^2;
   
    beta_1 = min(S_0,S_1);
%     if smoothness == 1
%     beta_2 =  (h/(6*c3))*((c2 + 2*c3*xph)^3 - (c2 + 2*c3*xmh)^3) + ((4*h^3))*c3^2*(xph - xmh); 
%     end

    if smoothness == 1
        beta_2 = (d'*A(:,:,i)*d);
    end
    
    if smoothness == 2
    beta_2 = (121*h/(600*c3))*((c2 + 2*c3*xph)^3 - (c2 + 2*c3*xmh)^3) + ((121*h^3)/5)*c3^2*(xph - xmh); 
    end
    

    %Note: Using smoothness indicator from q only
        
    tau_2 = abs(beta_2 - beta_1);

    omega_bar_12 = gam_12*(1 + tau_2/(epsilon + beta_1));
    omega_bar_22 = gam_22*(1 + tau_2/(epsilon + beta_2));

    omega_12 = omega_bar_12/(omega_bar_12 + omega_bar_22);
    omega_22 = omega_bar_22/(omega_bar_12 + omega_bar_22);
    wi = @(x) omega_12*p1(x) + omega_22*p2(x);

    uphm(i) = wi(xph);
    umhp(i) = wi(xmh);
end
end

function fhat = LF_flux(um,up,xmin)

    fm = (um).^2/2; fp = (up).^2/2;
    fhat = 0.5*( fm+fp -  xmin*(up-um) ); 
    return
end

function fhat = lf(a,uM,uP)
    fM = flux(uM); 
    fP = flux(uP);
    fhat = (1/2)*(fM+fP - a*(uP-uM));
end

function f = flux(u)
    f = (1/2)*u.^2;
end

function W = condensed_constrained_weight_matrix(k,m,r,i,x_array)
x_last_ghost = x_array(end) + x_array(2:k + 2);
x_first_ghost = x_array(1) - (x_array(end) - x_array(end - (k - 1):end-1));
x_ghost = [x_first_ghost; x_array; x_last_ghost];
dx_ghost = x_ghost(2:end) - x_ghost(1:end-1);

mod_n = @(i) i + k - 1;

%Build matrix for least squares
A = zeros(m,k);

for l = 1:k
    for j = 1:m
        A(j,l) = (x_ghost(mod_n(i - r + j))^l - x_ghost(mod_n(i - r + j - 1))^l)/(l * dx_ghost(mod_n(i - r + j - 1)));
    end
end


%Build L matrix
L = zeros(1,k);
for l = 1:k
    L(l) = (x_ghost(mod_n(i +1))^l - x_ghost(mod_n(i))^l)/(l * dx_ghost(mod_n(i)));
end

 %   %Can also add scaling factor for condition number
         %sc = 1/(m*max(max(abs(A))));

        %Solving for W with QR factorization
        [Q,R] = qr(A,0);
        W = R\(Q');          %(A'*A)\A'
        P = R \ (R' \ L');   %(A'*A)\(L');
        mu = -L*P;
        E = W + (P*L*W)/mu;
        F = P/mu;

        W = [E, F];
end

function A = stiffness_weight_matrix(k,m,r,i,x,W_i);
    dx = x(2:end) - x(1:end-1);
    h = dx(i);
    W = squeeze(W_i);
    W_2 = W(2,:); W_3 = W(3,:);
    xip1 = x(i + 1); xim1 = x(i);

    A_1 = h*(W_2'*W_2)*(xip1 - xim1);
    A_2 = 2*h*(W_2'*W_3)*(xip1^2 - xim1^2);
    A_3 = 4*h/3*(W_3' * W_3)*(xip1^3 - xim1^3);
    A_4 = 4*h^3*(W_3' * W_3)*(xip1 - xim1);

    A = A_1 + A_2 + A_3 + A_4;
end

function W = get_weight_matrix(k,m,r,i,x_array)
%Computes polynomial with degree on 
    %for k and i index with periodic boundaries
x_last_ghost = x_array(end) + x_array(2:k + 2);
x_first_ghost = x_array(1) - (x_array(end) - x_array(end - (k - 1):end-1));

x_ghost = [x_first_ghost; x_array; x_last_ghost];
dx_ghost = x_ghost(2:end) - x_ghost(1:end-1);

mod_n = @(i) i + k - 1;

%Build matrix for least squares
A = zeros(m,k);

for l = 1:k
    for j = 1:m
        A(j,l) = (x_ghost(mod_n(i - r + j))^l - x_ghost(mod_n(i - r + j - 1))^l)/(l * dx_ghost(mod_n(i - r + j - 1)));
    end
end

[Q,R] = qr(A);

W = R\(Q');
end

function order_plot(variable_name,L2Error,fig_num,Method,p,dx,err_name, k)
    figure(fig_num)
    hold on;
    plot(log10(dx),log10(L2Error),'.-','LineWidth',2,'MarkerSize',20)
    plot(log10(dx),p*log10(dx),'k','LineWidth',3) %Reference Slope = 3
    xlabel('Log_{10}(\Deltax)')
    ylabel('Log_{10}(Error)')
    title(sprintf("%s Convergence Plot using %s with k = %s", [err_name Method k]))
    legend('', sprintf('p=%i', p),'location','southeast')
end

function T = order_table(dx,n_array, L2Error, err_name)
    L2Order = (log10(L2Error(1:end-1))-log10(L2Error(2:end)))./(log10(dx(1:end-1))-log10(dx(2:end))) ;
    %Format table
    fractions = arrayfun(@(x) sprintf('1/%d', x), n_array, 'UniformOutput', false);
    T = table(fractions', L2Error(:),[0;L2Order(:)], ...
        'VariableNames',{'dx',sprintf('%s Error', err_name),sprintf('%s Order', err_name),});
    disp(T)
end


function int_res = one_dim_gauss_quad(fun,vals,weights)
      int_vals = fun(vals); %Evaluate interior function at vals
      int_res = sum(int_vals.*weights); %Obtain integral via  Legendre-Gauss quadrature
end




    