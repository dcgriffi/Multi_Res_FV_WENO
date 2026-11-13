%Multi-resolution WENO implementation with uniform rectangular mesh (constrained) 
    %%% Run pde = linear for order or pde = burgers for behavior around shock

clear all; close all;

pth = fullfile(pwd, 'Results');  % figures/tables are saved to ./Results/2D
if ~exist(fullfile(pth, '2D'), 'dir'), mkdir(fullfile(pth, '2D')); end
%pth = '';

%Uncomment the desired time
%t_final = 0.3; %Before shock 1D
%t_final =2; %After shock 1D
%t_final = 0.3/2; %Before shock 2D
t_final = 2/2; %After shock 2D

avoid_shock = true;

if  t_final == 2 || t_final == (2/2) 
    avoid_shock = true;
end

if t_final == 0.3 || t_final == (0.3/2)
    avoid_shock = false;
end

p = 3;

pde = "burgers";
%pde = "linear";
alpha = 1/3; beta = 2/3; 
u_initial = @(x,y) alpha + beta*sin((x+y));
%u_initial = @(x,y) alpha + beta*sin((x));

if pde == "burgers"
f = @(u) u.^2/2;
g = @(u) u.^2/2;
elseif pde == "linear"
f = @(u) u;
g = @(u) u; 
end

n_array=[20,40,80,160];

%Preallocate arrays for memory
DX = zeros(length(n_array),1); dt = zeros(length(n_array),1); 
L2_Error_quadratic = zeros(length(n_array),1); L1Error = zeros(length(n_array),1); 
L2_Error_linear = zeros(length(n_array),1);
LinfError= zeros(length(n_array),1); err = zeros(length(n_array),1); 

%xL = 0; xR = 4*pi;
xL = 0; xR = 2*pi;

yL = xL; yR = xR;
   %Cycle through all ranges of dx
   for kk=1:length(n_array)

        %Build spatial array
        N = n_array(kk);
        DX(kk)=(xR-xL)/N;
        h=DX(kk); 
        x = (xL:h:xR)';
        y = (yL:h:yR); 
        dx = x(2:end) - x(1:end-1);
        dy = y(2:end) - y(1:end-1);

        nX = length(x)-1;
        nY = length(y)-1;
        N = (nX)*(nY);
       
        [ind_array,sten_array, coord_array, sten_coord_array] = get_rectangle_mesh_info(nX,nY,x,y,dx,dy,N);
        v_vec = build_cell_avg_vec(coord_array,u_initial,N);

        %Precompute quad_vals and length of sides
        [quad_vals, len_vec] = get_quad_info(coord_array,N);

      
        %Third order approximation
        k = 6; m = 9;

        %Precompute weight matrix for coefficents
        W = NaN(k, m+1, N);
        for i = 1:N
            W(:,:,i) = rectangle_condensed_constrained_weight_matrix(k,m,i,coord_array,sten_coord_array);
        end

        B = NaN(m+1, m+1, N);
        for i = 1:N
            B(:,:,i) = rectangle_smoothness_matrix(squeeze(W(:,:,i)),i,coord_array);
        end



        %Change both functions into plottable array
        t = 0;
        while t < t_final     
            %Choose appropriate time step
             h = max(dx);
             a = max(abs(v_vec));
             dt = (h)/(2*a); 
             %dt = h*0.07;

            if t_final - t < dt
                 dt = t_final - t;
            end
    
            %RK3

            u_1 = v_vec + dt*L(v_vec,W,B, quad_vals,len_vec,sten_array,coord_array,sten_coord_array,N,f,g,a);
            u_2 = 3/4*v_vec +1/4*u_1 + 1/4*dt*L(u_1,W,B, quad_vals,len_vec,sten_array,coord_array,sten_coord_array,N,f,g,a);
            u_next = 1/3*v_vec + 2/3*u_2 + 2/3*dt*L(u_2,W,B, quad_vals,len_vec,sten_array,coord_array,sten_coord_array,N,f,g,a);
    
            v_vec = u_next;
            t = t + dt;
        end
        
        %Prep for visualization

        if pde == "burgers"
            u_ref= @(x,y) solveBurgers2(alpha, beta, 2*t_final,(x+y));
            %u_ref= @(x,y) solveBurgers2(alpha, beta, t_final,(x));
            %u_ref= @(x,y) solveBurgers2(alpha, beta, t_final,(y));


         elseif pde == "linear"
            u_ref = @(x,y) u_initial(x - t_final, y - t_final);

        end

        variable_name = "1/3 + 2/3sin(x+y)"; 
        Method = "Least Squares Poly Approx"; 

        

        v_array = reshape(v_vec,nX,nY);
        u_exact_avg = build_cell_avg_vec(coord_array,u_ref,N);
        u_exact_array = reshape(u_exact_avg, nX,nY);
        x_avg = (x(1:nX) + x(2:nX+1))/2;
        y_avg = (y(1:nY) + y(2:nY+1))/2;
        
        %Publication-style profile along the diagonal (N = 80): all cell
        %averages plotted against x + y (mod 2*pi), with the exact solution
        %as a red line. Saves a vector PDF to Results/2D.
        if nX == 80 && pde == "burgers"
            plot_diagonal_profile(x_avg, y_avg, v_array, ...
                @(s) solveBurgers2(alpha, beta, 2*t_final, s), ...
                fullfile(pth, '2D', sprintf('rect_burgers_t%s_N%d.pdf', strrep(num2str(t_final), '.', 'p'), nX)));
        end
        %Errors in the cell averages, area-weighted and normalized:
        %  L1 = sum |V_j||e_j| / sum |V_j|,  L2 = (sum |V_j| e_j^2 / sum |V_j|)^(1/2),
        %  Linf = max |e_j|,
        %over all cells, or over the cells away from the shock when avoid_shock.
        uDif = v_vec - u_exact_avg;
        cell_area = (coord_array(:,2) - coord_array(:,1)).*(coord_array(:,4) - coord_array(:,3));
        if avoid_shock
            mask = away_from_shock(N, coord_array, t_final, alpha);
        else
            mask = true(N,1);
        end
        w = cell_area(mask)/sum(cell_area(mask));
        L1Error(kk) = sum(w.*abs(uDif(mask)));
        L2_Error_quadratic(kk) = sqrt(sum(w.*uDif(mask).^2));
        LinfError(kk) = max(abs(uDif(mask)));

   end

    format shorte
    p = 3;
    fig_num = 5; err_name = "L2 Error";
    T1 = order_table(DX,n_array, L1Error, "L1 Error");
    T2 = order_table(DX,n_array, L2_Error_quadratic, "L2 Error");
    T3 = order_table(DX,n_array, LinfError, "Linf Error");



    
%%%%%%Start of functions for multires%%%%%%%%%%%%
function u_hold = L(v_vec,W,B, quad_vals,len_vec,sten_array,coord_array,sten_coord_array,N,f,g,a)
    quad_res = weno_reconstruction(W,B,sten_array,coord_array,sten_coord_array,quad_vals,v_vec,N); 
    flux_array = find_fluxes(quad_res,sten_array,N,f,g,a);

    u_hold = NaN(N,1);
    for i = 1:N
        x_m = coord_array(i,1); x_p = coord_array(i,2); 
        y_m = coord_array(i,3); y_p = coord_array(i,4);
        sum = 0;
        for e = 1:4
            weights = [0.5,0.5];
            len = len_vec(i,e);
            sum = sum + len*(weights(1)*flux_array(i,e,1) + weights(2)*flux_array(i,e,2));
        end
        area = (x_p - x_m)*(y_p - y_m);
        u_hold(i) = -(1/area)*sum;
    end
end

function quad_res = weno_reconstruction(W,B,sten_array,coord_array,sten_coord_array,quad_vals,v_vec,N)
   quad_res = NaN(N,4,2);
   gam_12 = 1/11; gam_22 = 10/11;  epsilon = 1e-5;
   m=9;
   for i = 1:N
        stencil = sten_array(i,:);
        d = zeros(m+1,1);
        for j = 1:m
            cell = stencil(j);
            d(j) = v_vec(cell);
        end
        d(m+1) = -v_vec(i);
    
        C = W(:,:,i)*d;


        q1 = @(x,y) v_vec(i);

        q2 = @(x,y) C(1) + C(2)*x + C(3)*y+ C(4)*x.*y + C(5)*x.^2 + C(6)*y.^2;
    
        p1 = @(x,y) q1(x,y);
        p2 = @(x,y) (11/10)*q2(x,y) - (1/10)*p1(x,y);
    
        %Build B_1 (Pick minimum from four different stencils)

        x_m = sten_coord_array(i,1,1); x_p = sten_coord_array(i,1,2);
        y_m = sten_coord_array(i,1,3); y_p = sten_coord_array(i,1,4);

        area = (x_p - x_m)*(y_p - y_m);

        stencil_1 = [1,2,4];
        stencil_2 = [1,4,6];
        stencil_3 = [1,6,8];
        stencil_4 = [1,8,2];

        [C2_1, C3_1, b1_1] = compute_lin_smooth(i,sten_coord_array,sten_array,v_vec, stencil_1);
        [C2_2, C3_2, b1_2] = compute_lin_smooth(i,sten_coord_array,sten_array,v_vec, stencil_2);
        [C2_3, C3_3, b1_3] = compute_lin_smooth(i,sten_coord_array,sten_array,v_vec, stencil_3);
        [C2_4, C3_4, b1_4] = compute_lin_smooth(i,sten_coord_array,sten_array,v_vec, stencil_4);

        if b1_3 >= b1_1
            gam_bar_13 = 1;
        else 
            gam_bar_13 = 10;
        end

        gam_bar_11 = 11 - gam_bar_13;

        gam_bar_12 = (gam_bar_11 + gam_bar_13)/2;
        gam_bar_14 = gam_bar_12;

        bar_sum = (gam_bar_11 + gam_bar_12 + gam_bar_13 + gam_bar_14);
        gam_11 = gam_bar_11/(bar_sum); gam_12 = gam_bar_12/(bar_sum);
        gam_13 = gam_bar_13/(bar_sum); gam_14 = gam_bar_14/(bar_sum);


        beta_dif_avg =(abs(b1_1 - b1_2)+abs(b1_1 - b1_3) + abs(b1_1 - b1_4) + ...
                     + abs(b1_2 - b1_3)+abs(b1_2 - b1_4) + abs(b1_3 - b1_4))/6;
        % %1D approach 
        % sigma_1 = (gam_11)*(1 + (beta_dif_avg)/(b1_1 + epsilon));
        % sigma_2 = (gam_12)*(1 + (beta_dif_avg)/(b1_2 + epsilon));
        % sigma_3 = (gam_13)*(1 + (beta_dif_avg)/(b1_3 + epsilon));
        % sigma_4 = (gam_14)*(1 + (beta_dif_avg)/(b1_4 + epsilon));
        % 
        % sigma = sigma_1 + sigma_2  + sigma_3 + sigma_4;
        % 
        % 
        % beta_1 = (1/sigma^2)*(sigma_1*(b1_1)^(1/2) + sigma_2*b1_2^(1/2)+...
        %                     + sigma_3*b1_3^(1/2)   + sigma_4*b1_4^(1/2))^2;

        %2D Lit approach
        sigma_1 = 1/4*(1 + (beta_dif_avg)^2/(b1_1 + epsilon));
        sigma_2 = 1/4*(1 + (beta_dif_avg)^2/(b1_2 + epsilon));
        sigma_3 = 1/4*(1 + (beta_dif_avg)^2/(b1_3 + epsilon));
        sigma_4 = 1/4*(1 + (beta_dif_avg)^2/(b1_4 + epsilon));

        sigma = sigma_1 + sigma_2  + sigma_3 + sigma_4;


        beta_1 = (area/sigma^2)*(sigma_1*(C2_1) + sigma_2*(C2_2)+...
                               + sigma_3*(C2_3) + sigma_4*(C2_4))^2+...
                +(area/sigma^2)*(sigma_1*(C3_1) + sigma_2*(C3_2)+...
                               + sigma_3*(C3_3) + sigma_4*(C3_4))^2;
        

        %2D updated approach
        % sigma_1 = (gam_11)*(1 + (beta_dif_avg)^2/(b1_1 + epsilon));
        % sigma_2 = (gam_12)*(1 + (beta_dif_avg)^2/(b1_2 + epsilon));
        % sigma_3 = (gam_13)*(1 + (beta_dif_avg)^2/(b1_3 + epsilon));
        % sigma_4 = (gam_14)*(1 + (beta_dif_avg)^2/(b1_4 + epsilon));
        % 
        % 
        % sigma = sigma_1 + sigma_2  + sigma_3 + sigma_4;
        % 
        % 
        % beta_1 = (area/sigma^2)*(sigma_1*(C2_1) + sigma_2*(C2_2)+...
        %                        + sigma_3*(C2_3) + sigma_4*(C2_4))^2+...
        %         +(area/sigma^2)*(sigma_1*(C3_1) + sigma_2*(C3_2)+...
        %                        + sigma_3*(C3_3) + sigma_4*(C3_4))^2;
        % 

       
        %Minimum approach
        %beta_1 = min([b1_1,b1_2,b1_3,b1_4]);
        %beta_1 = 0;

        gam_12 = 1/11; gam_22 = 10/11;  %epsilon = 1e-5;

        %Build B_2
        beta_2 = d'*B(:,:,i)*d;
    
        %Note: Using smoothness indicator from q only

        tau_2 = abs(beta_2 - beta_1);
    
        omega_bar_12 = gam_12*(1 + tau_2/(epsilon + beta_1));
        omega_bar_22 = gam_22*(1 + tau_2/(epsilon + beta_2));
    
        omega_12 = omega_bar_12/(omega_bar_12 + omega_bar_22);
        omega_22 = omega_bar_22/(omega_bar_12 + omega_bar_22);
        
        u_cell = @(x,y) omega_12*p1(x,y) + omega_22*p2(x,y);

        x_m = coord_array(i,1); x_p = coord_array(i,2); 
        y_m = coord_array(i,3); y_p = coord_array(i,4);
        
        %Edge 1
        quad_res(i,1,:) = u_cell(x_p, quad_vals(i,1,:));

        %Edge 2
        quad_res(i,2,:) = u_cell(quad_vals(i,2,:),y_p);

        %Edge 3
        quad_res(i,3,:) = u_cell(x_m,quad_vals(i,3,:));

        %Edge 4
        quad_res(i,4,:) = u_cell(quad_vals(i,4,:),y_m);
    end
end

function [C2, C3, b1] = compute_lin_smooth(i,sten_coord_array,sten_array,v_vec, stencil)

        x_m = sten_coord_array(i,1,1); x_p = sten_coord_array(i,1,2);
        y_m = sten_coord_array(i,1,3); y_p = sten_coord_array(i,1,4);

        area = (x_p - x_m)*(y_p - y_m);

        get_cen = @(x1,x2) (x1 + x2)/2;

        xc_1 = get_cen(x_m,x_p); yc_1 = get_cen(y_m,y_p);

        v1 = v_vec(i);

        x_m = sten_coord_array(i,stencil(2),1); x_p = sten_coord_array(i,stencil(2),2);
        y_m = sten_coord_array(i,stencil(2),3); y_p = sten_coord_array(i,stencil(2),4);

        xc_2 = get_cen(x_m,x_p); yc_2 = get_cen(y_m,y_p);

        cell_2 = sten_array(i,stencil(2));
        v2 = v_vec(cell_2);

        x_m = sten_coord_array(i,stencil(3),1); x_p = sten_coord_array(i,stencil(3),2);
        y_m = sten_coord_array(i,stencil(3),3); y_p = sten_coord_array(i,stencil(3),4);

        xc_3 = get_cen(x_m,x_p); yc_3 = get_cen(y_m,y_p);

        cell_3 = sten_array(i,stencil(3));
        v3 = v_vec(cell_3);

        dx2 = xc_2 - xc_1 ; dx3 = xc_3 - xc_1;
        dy2 = yc_2 - yc_1;  dy3 = yc_3 - yc_1;
        dv2 = v2 - v1;      dv3 =  v3 - v1;

        C3 = (dx2*dv3 - dx3*dv2)/(dy3*dx2 - dx3*dy2);

        if dx3 == 0
            C2 = (dv2 - C3*dy2)/dx2;
        elseif dx2 == 0
            C2 = (dv3 - C3*dy3)/dx3;
        else
            C2 = (dy3*dv2 - dv3*dy2)/(dx2*dy3 - dy3*dx3);
        end

        b1 = (C2^2+C3^2)*area;
end

function B = rectangle_smoothness_matrix(W,i,coord_array)

    x_m = coord_array(i,1); x_p = coord_array(i,2); y_m = coord_array(i,3); y_p = coord_array(i,4);
    area = (x_p - x_m)*(y_p - y_m);

    L1 = (x_p - x_m)*(y_p - y_m);
    L2 = (x_p^2 - x_m^2)*(y_p-y_m)/2;
    L3 = (y_p^2 - y_m^2)*(x_p-x_m)/2;
    L4 = (x_p^2 - x_m^2)*(y_p^2 - y_m^2)/4;
    L5 = (x_p^3 - x_m^3)*(y_p-y_m)/3;
    L6 = (y_p^3 - y_m^3)*(x_p-x_m)/3;

    W2 = W(2,:); W3 = W(3,:); 
    W4 = W(4,:); W5 = W(5,:); W6 = W(6,:);

    B1 = ((W2')*W2 + (W3')*W3 + area*(W4'*W4 + 4*(W5')*W5 + 4*(W6')*W6))*L1;
    B2 = (4*(W2')*W5 + 2*(W3')*W4)*L2;
    B3 = (2*(W2')*W4 + 4*(W3')*W6)*L3;
    B4 = (4*(W4')*W5 + 4*(W4')*W6)*L4;
    B5 = (4*(W5')*W5 + (W4'*W4))*L5;
    B6 = ((W4')*W4 + 4*(W6')*W6)*L6;

   % B = 121/100*(B1 + B2 + B3 + B4 + B5 + B6) ;
    B = (B1 + B2 + B3 + B4 + B5 + B6);

end


function flux_array = find_fluxes(quad_res,sten_array,N,f,g,a)
    
    flux_array = NaN(N,4,2);
    
    %Find fluxes for edge 1 in each cell
    for i = 1:N
        %Obtain u_m from current cell (for edge 1)
        u_m = quad_res(i,1,:);

        %Obtain u_p from adjacent cell (2) correct edge (3)
        adj_cell = sten_array(i,2);
        u_p = quad_res(adj_cell, 3,:);

        %Compute flux with norm array pointed out of current cell
        norm = [1,0];
        flux_array(i,1,1) = lf(u_m(1),u_p(1),norm,f,g,a);
        flux_array(i,1,2) = lf(u_m(2),u_p(2),norm,f,g,a);
    end

    %Find fluxes for edge 2 in each cell
    for i = 1:N
        %Obtain u_m from current cell (for edge 2, horizontal)
        u_m = quad_res(i,2,:);

        %Obtain u_p from adjacent cell (4) correct edge (4)
        adj_cell = sten_array(i,4);
        u_p = quad_res(adj_cell, 4,:);

        %Compute flux with norm array pointed out of current cell
        norm = [0,1];
        flux_array(i,2,1) = lf(u_m(1),u_p(1),norm,f,g,a);
        flux_array(i,2,2) = lf(u_m(2),u_p(2),norm,f,g,a);
    end

    %Assign fluxes to correct edges (negative copies to 3 and 4)
    for i = 1:N
        adj_cell = sten_array(i,6);
        flux_array(i,3,:) = -flux_array(adj_cell,1,:);

        adj_cell = sten_array(i,8);
        flux_array(i,4,:) = -flux_array(adj_cell,2,:);
    end

end

function fhat = lf(u_m,u_p,norm,f,g,a)
    F =    NaN(2,1);
    f_m = f(u_m); f_p = f(u_p);
    F(1) = (f_m+f_p);
    g_m = g(u_m); g_p = g(u_p);
    F(2) = (g_m+g_p);
    fhat = (1/2)*(dot(F,norm) - a*(u_p-u_m));
end   
    
%%%%%%Functions from Least Square
function T = order_table(dx,n_array, L2Error, err_name)
    L2Order = (log10(L2Error(1:end-1))-log10(L2Error(2:end)))./(log10(dx(1:end-1))-log10(dx(2:end))) ;
    %Format table
    fractions = arrayfun(@(x) sprintf('1/%d', x), n_array, 'UniformOutput', false);
    T = table(fractions', L2Error(:),[0;L2Order(:)], ...
        'VariableNames',{'dx',sprintf('%s Error', err_name),sprintf('%s Order', err_name),});
    disp(T)
end


function mask = away_from_shock(N, coord_array, t_final, alpha)
    %Cells whose x + y range lies at distance more than 1 from the shock.
    %The 2D solution is w(x + y, 2t), so the shock is at x + y = pi + 2*alpha*t
    %(and its periodic copy at x + y = 3*pi + 2*alpha*t).
    xShock = 2*pi + alpha * 2*t_final;
    diag_shock1 = xShock - pi;
    diag_shock2 = xShock + pi;
    low_1 = diag_shock1 - 1;  high_1 = diag_shock1 + 1;
    low_2 = diag_shock2 - 1;  high_2 = diag_shock2 + 1;
    mask = false(N,1);
    for j = 1:N
        x_m = coord_array(j,1); x_p = coord_array(j,2); y_m = coord_array(j,3); y_p = coord_array(j,4);
        mask(j) = ((x_p + y_p) <= low_1) || (((x_m + y_m) >= high_1) && ((x_p + y_p) <= low_2)) || ((x_m + y_m) >= high_2);
    end
end


function result = Gauss_quad_2D(xlow,xhigh,ylow,yhigh,interior)
        [x_vals,x_weights] = lgwt(2,xlow, xhigh);
        [y_vals, y_weights] = lgwt(2,ylow, yhigh);

        %integrand = interior(x_vals,y_vals');
        
        integrand = NaN(length(x_vals),length(y_vals));
        for p=1:length(x_vals)
            for q=1:length(y_vals)
                integrand(p,q) = interior(x_vals(p), y_vals(q));
            end
        end


        result =  sum(sum(integrand.*x_weights.*y_weights'));
end

function W = rectangle_condensed_constrained_weight_matrix(k,m,i,coord_array,sten_coord_array)
        %Center the least squares on cell i and scale by sqrt(area). Building
        %the system in absolute coordinates makes A'A badly conditioned.
        xc = (coord_array(i,1) + coord_array(i,2))/2; yc = (coord_array(i,3) + coord_array(i,4))/2;
        s = sqrt((coord_array(i,2) - coord_array(i,1))*(coord_array(i,4) - coord_array(i,3)));


         %Build matrix for least squares
        A = zeros(m,k);
        for j = 1:m
            x_m = (sten_coord_array(i,j,1)-xc)/s; x_p = (sten_coord_array(i,j,2)-xc)/s; y_m = (sten_coord_array(i,j,3)-yc)/s; y_p = (sten_coord_array(i,j,4)-yc)/s;
            area = (x_p - x_m)*(y_p - y_m);
            A(j,1) = 1/area*(x_p - x_m)*(y_p - y_m);
            A(j,2) = 1/area*(x_p^2 - x_m^2)*(y_p-y_m)/2;
            A(j,3) = 1/area*(y_p^2 - y_m^2)*(x_p-x_m)/2;
            A(j,4) = 1/area*(x_p^2 - x_m^2)*(y_p^2 - y_m^2)/4;
            A(j,5) = 1/area*(x_p^3 - x_m^3)*(y_p-y_m)/3;
            A(j,6) = 1/area*(y_p^3 - y_m^3)*(x_p-x_m)/3;
        end
        
        %Build L matrix
        L = zeros(1,k);
        cell = i;
        x_m = (coord_array(cell,1)-xc)/s; x_p = (coord_array(cell,2)-xc)/s; y_m = (coord_array(cell,3)-yc)/s; y_p = (coord_array(cell,4)-yc)/s;
        area = (x_p - x_m)*(y_p - y_m);
        L(1) = 1/area*(x_p - x_m)*(y_p - y_m);
        L(2) = 1/area*(x_p^2 - x_m^2)*(y_p-y_m)/2;
        L(3) = 1/area*(y_p^2 - y_m^2)*(x_p-x_m)/2;
        L(4) = 1/area*(x_p^2 - x_m^2)*(y_p^2 - y_m^2)/4;
        L(5) = 1/area*(x_p^3 - x_m^3)*(y_p-y_m)/3;
        L(6) = 1/area*(y_p^3 - y_m^3)*(x_p-x_m)/3;
        
          
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

        %Map coefficients of 1, z, w, zw, z^2, w^2 (z = (x-xc)/s, w = (y-yc)/s)
        %back to coefficients of 1, x, y, xy, x^2, y^2, so the rest of the
        %code is unchanged.
        W = scaled_to_absolute_2d(xc,yc,s)*W;





end


  function [ind_array,sten_array, coord_array, sten_coord_array] = get_rectangle_mesh_info(nX,nY,x,y,dx,dy,N)
 
        %Build array to enumerate cells (increasing left to right, then top to
        %bottom)
        
        ind = 0;
        ind_array = zeros(nX,nY);
        for n = 1:nY
            for j = 1:nX
                ind = ind + 1;
                ind_array(j,n) = ind;
            end
        end
      
        %Pad array with ghost cells for periodicity
        pad_ind_array = [ind_array(:,end), ind_array, ind_array(:,1) ];
        pad_ind_array = [pad_ind_array(end,:); pad_ind_array; pad_ind_array(1,:) ];

        pad_x = [x(1) - dx(end); x; x(end) + dx(1)];
        pad_y = [y(1) - dy(end), y, y(end) + dy(1)];
        
        %Build array to store stencil for each i
        ind = 0;
        sten_array = NaN(N,9);
        for yy = 2:nY+1
            for xx = 2:nX+1
                ind = ind + 1;
                sten_array(ind,1) = pad_ind_array(xx,  yy);
                sten_array(ind,2) = pad_ind_array(xx+1,yy);
                sten_array(ind,3) = pad_ind_array(xx+1,yy+1);
                sten_array(ind,4) = pad_ind_array(xx,  yy+1);
                sten_array(ind,5) = pad_ind_array(xx-1,yy+1);
                sten_array(ind,6) = pad_ind_array(xx-1,yy);
                sten_array(ind,7) = pad_ind_array(xx-1,yy-1);
                sten_array(ind,8) = pad_ind_array(xx,  yy-1);
                sten_array(ind,9) = pad_ind_array(xx+1,yy-1);
            end
        end


                
        ind = 0;
        coord_array = NaN(N,4);
        for n = 1:nY
            for j = 1:nX
                ind = ind + 1;
                coord_array(ind,1) = x(j);
                coord_array(ind,2) = x(j+1);
                coord_array(ind,3) = y(n);
                coord_array(ind,4) = y(n+1);
            end
        end

        sten_coord_array = NaN(N,9,4);
        ind = 0;
         for yy = 3:nY+2
            for xx = 3:nX+2
                ind = ind + 1;
                sten_coord_array(ind,1,:) = [pad_x(xx-1),pad_x(xx),  pad_y(yy-1),pad_y(yy)];
                sten_coord_array(ind,2,:) = [pad_x(xx),  pad_x(xx+1),pad_y(yy-1),pad_y(yy)];
                sten_coord_array(ind,3,:) = [pad_x(xx),  pad_x(xx+1),pad_y(yy),  pad_y(yy+1)];
                sten_coord_array(ind,4,:) = [pad_x(xx-1),pad_x(xx),  pad_y(yy),  pad_y(yy+1)];
                sten_coord_array(ind,5,:) = [pad_x(xx-2),pad_x(xx-1),pad_y(yy),  pad_y(yy+1)];
                sten_coord_array(ind,6,:) = [pad_x(xx-2),pad_x(xx-1),pad_y(yy-1),pad_y(yy)];
                sten_coord_array(ind,7,:) = [pad_x(xx-2),pad_x(xx-1),pad_y(yy-2),pad_y(yy-1)];
                sten_coord_array(ind,8,:) = [pad_x(xx-1),pad_x(xx),  pad_y(yy-2),pad_y(yy-1)];
                sten_coord_array(ind,9,:) = [pad_x(xx),  pad_x(xx+1),pad_y(yy-2),pad_y(yy-1)];
            end
         end
  end


  function [quad_vals, len_vec] = get_quad_info(coord_array,N)

        quad_vals = NaN(N,4,2);
        len_vec = NaN(N,4);
        for i = 1:N
          x_m = coord_array(i,1);  x_p = coord_array(i,2); 
          y_m = coord_array(i,3);  y_p = coord_array(i,4);
          
          dy = y_p - y_m;
          dx = x_p - x_m;

          %Edge 1 & 3
          [vals,weights] = lgwt(2,y_m,y_p); 
          quad_vals(i,1,:) = vals;

          quad_vals(i,3,:) = vals;

          len_vec(i,1) = dy;
          len_vec(i,3) = dy;

          %Edge 2 & 4
          [vals,weights] = lgwt(2,x_m,x_p); 
          quad_vals(i,2,:) = vals;

          quad_vals(i,4,:) = vals;

          len_vec(i,2) = dx;
          len_vec(i,4) = dx;  
        end
  end

 function v_vec = build_cell_avg_vec(coord_array,u_initial,N)

    v_vec = zeros(N,1);
    %Initialize cell average and fill in initial conditions
    for j=1:N
        x_m = coord_array(j,1); x_p = coord_array(j,2); y_m = coord_array(j,3); y_p = coord_array(j,4);
        area = (x_p - x_m)*(y_p - y_m);
        v_vec(j) = Gauss_quad_2D(x_m,x_p,y_m,y_p,u_initial)/area;
    end 
 end


 function plot_2D_rectangle_exact_approx(x,y,nX,nY,C,ind_array,u_initial,fig_num,h,pth)
    figure(fig_num)
    x_dense = NaN(30,1);
    y_dense = NaN(1,30);
    u_approx_array = NaN(30,30);
    u_exact_array = NaN(30,30);
    int = 1;
    for i = 1:nY
        int_points_y = (y(i):h/2:y(i+1));
         for m = 1:3
             y_point = int_points_y(m);
             for j = 1:nX           
                 cell = ind_array(j,i);
                 u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
                 int_points_x = (x(j):h/2:x(j+1)); 
                 x_dense((1+3*(j-1)):(3+3*(j-1))) = int_points_x;
                 y_dense(m+3*(i - 1)) = y_point;
                 u_approx_array((1+3*(j-1)):(3+3*(j-1)),(m+3*(i - 1))) = u_cell(int_points_x,y_point);
                 u_exact_array((1+3*(j-1)):(3+3*(j-1)),(m+3*(i - 1))) = u_initial(int_points_x,y_point);
                 
             end
         end
    end
    [X,Y] = meshgrid(x_dense,y_dense);
    surf(X,Y,u_approx_array);
    title(sprintf("2D Plot of multiresolution results for 1/3 + 2/3 sin(x + y)"))
    
    fig_num = fig_num + 1;
    figure(fig_num)
    surf(X,Y,u_exact_array);
    title(sprintf("2D Plot of exact solution for 1/3 + 2/3 sin(x + y)"))
 end


 function difference = test_cell_avg_difference_rectangle(coord_array,nX,nY,N,C,v_array);
 %Check the difference between the function and cell averages for
        %each cell (it works)
        
        cell_avg = zeros(N,1);
        for i = 1:N
            cell = i;
            %Compute cell average using polynomial c
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
            x_m = coord_array(i,1); x_p = coord_array(i,2); y_m = coord_array(i,3); y_p = coord_array(i,4);
            area = (x_p - x_m)*(y_p- y_m);
            %Obtain local error estimate with quadrature
            interior= @(x,y) u_cell(x,y);
            cell_avg(i) = Gauss_quad_2D(x_m,x_p,y_m,y_p,interior)/area;

        end
        

        cell_avg_array = reshape(cell_avg, [nX,nY]);
        difference = abs(v_array - cell_avg_array);
 end

 function int_res = one_dim_gauss_quad(fun,vals,weights)
      int_vals = fun(vals); %Evaluate interior function at vals
      int_res = sum(int_vals.*weights); %Obtain integral via  Legendre-Gauss quadrature
 end

function T = scaled_to_absolute_2d(xc,yc,s)
    %Rows: coefficients of 1, x, y, xy, x^2, y^2. Columns: coefficients of
    %1, z, w, zw, z^2, w^2 with z = (x-xc)/s and w = (y-yc)/s.
    T = [1, -xc/s, -yc/s, xc*yc/s^2,  xc^2/s^2,  yc^2/s^2;
         0,  1/s,   0,    -yc/s^2,   -2*xc/s^2,  0;
         0,  0,     1/s,  -xc/s^2,    0,        -2*yc/s^2;
         0,  0,     0,     1/s^2,     0,         0;
         0,  0,     0,     0,         1/s^2,     0;
         0,  0,     0,     0,         0,         1/s^2];
end


function plot_diagonal_profile(x_avg, y_avg, v_array, w_exact, fname)
    %Cell averages against s = x + y (mod 2*pi). For data that depend only on
    %x + y, all cells collapse onto one curve, and any spread shows
    %oscillations. Exact solution as a red line, cell averages as green
    %symbols, no title (details go in the LaTeX caption).
    [X, Y] = ndgrid(x_avg, y_avg);
    s = mod(X(:) + Y(:), 2*pi);
    sf = linspace(0, 2*pi, 2000);
    fig = figure('Units','inches','Position',[1 1 5.5 4.0],'Color','w');
    plot(sf, arrayfun(w_exact, sf), 'r-', 'LineWidth', 1.5); hold on
    plot(s, v_array(:), 'o', 'Color', [0 0.5 0], 'MarkerSize', 3, 'LineWidth', 0.75);
    xlim([0 2*pi]); ylim([-0.45 1.1]);
    xticks(0:pi/2:2*pi); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
    xlabel('x + y'); ylabel('u');
    legend({'Exact','Multi-resolution WENO'}, 'Location', 'southwest', 'Box', 'off');
    set(gca, 'FontSize', 13, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'on');
    exportgraphics(fig, fname, 'ContentType', 'vector');
end
