%FV_LS_2D_Rectangle.m
%   Finite volume scheme on a uniform rectangular grid with the constrained
%   quadratic least squares reconstruction only (no WENO weights), for the
%   2D Burgers' equation with u0 = 1/3 + 2/3 sin(x + y) on [0,2*pi]^2.
%   The 2D solution is w(x + y, 2t), so t = 0.15 and t = 1 correspond to
%   t = 0.3 and t = 2 in 1D; the shock forms at t = 0.75. Before the shock the
%   scheme is third order; after the shock the reconstruction oscillates.
%   Prints L1, L2 and Linf errors of the cell averages (area-weighted and
%   normalized; away from the shock after it forms) and saves a profile of
%   all cell averages against x + y for N = 80 to Results/2D.
%   Run pde = "linear" for an advection test.

clear; close all;
pth = fullfile(pwd, 'Results');  % figures are saved to ./Results/2D
if ~exist(fullfile(pth, '2D'), 'dir'), mkdir(fullfile(pth, '2D')); end
t_final = 0.3/2;        % 0.3/2 (before the shock) or 2/2 (after the shock)
avoid_shock = (t_final > 0.75);   % measure the error away from the shock after it forms

p = 3;

pde = "burgers";        % "burgers" or "linear"
alpha = 1/3; beta = 2/3; 
%u_initial = @(x,y) alpha + beta*sin((x+y)/2); 
u_initial = @(x,y) alpha + beta*sin((x+y)); 

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

% xL = 0; xR = 4*pi;
xL = 0; xR = 2*pi;

yL = xL; yR = xR;
   %Cycle through all ranges of dx
   for kk=1:length(n_array)

        %Build spatial array
        N = n_array(kk);
        DX(kk)=(xR-xL)/N;
        h=DX(kk); 
        x = (xL:h:xR)'; %We have nSteps + 1 points here
        y = (yL:h:yR); %We have nSteps + 1 points here
        dx = x(2:end) - x(1:end-1);
        dy = y(2:end) - y(1:end-1);

        nX = length(x)-1;
        nY = length(y)-1;
        N = (nX)*(nY);
       
        [ind_array,sten_array, coord_array, sten_coord_array] = get_rectangle_mesh_info(nX,nY,x,y,dx,dy,N);
        v_vec = build_cell_avg_vec(coord_array,u_initial,N);
%         v_array = reshape(v_vec, [nX,nY]);

        %Precompute quad_vals and quad_weights
        [quad_vals, quad_weights, len_vec] = get_quad_info(coord_array,N);

      
        %Third order approximation
        k = 6; m = 9;

        %Precompute weight matrix for coefficents
        W = zeros(k, m+1, N);
        for i = 1:N
            W(:,:,i) = rectangle_condensed_constrained_weight_matrix(k,m,i,coord_array,sten_coord_array);
        end

        %Change both functions into plottable array
        t = 0;
        while t < t_final     
            %Choose appropriate time step
             h = max(dx);
             a = max(abs(v_vec));
             dt = (h)/(4*a); 
             %dt = h*0.07;

            if t_final - t < dt
                 dt = t_final - t;
            end
    
            %RK3
            u_1 = v_vec + dt*L(v_vec,W, quad_vals,quad_weights,len_vec,sten_array,coord_array,N,f,g,a);
            u_2 = 3/4*v_vec +1/4*u_1 + 1/4*dt*L(u_1,W, quad_vals,quad_weights,len_vec,sten_array,coord_array,N,f,g,a);
            u_next = 1/3*v_vec + 2/3*u_2 + 2/3*dt*L(u_2,W, quad_vals,quad_weights,len_vec,sten_array,coord_array,N,f,g,a);
    
            v_vec = u_next;
            t = t + dt;
        end
        
        %Prep for visualization

        if pde == "burgers"
            %u_ref= @(x,y) solveBurgers2(alpha, beta, t_final,(x+y)./2);
            u_ref= @(x,y) solveBurgers2(alpha, beta, 2*t_final,(x+y));

         elseif pde == "linear"
            u_ref = @(x,y) u_initial(x - t_final, y - t_final);
            %u_ref = @(x,y) u_initial(x - t_final, y - t_final);

        end

        %Errors in the cell averages, area-weighted and normalized, over all
        %cells or over the cells away from the shock
        u_exact_avg = build_cell_avg_vec(coord_array,u_ref,N);
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

        %Profile of all cell averages against x + y (N = 80)
        if nX == 80 && pde == "burgers"
            x_avg = (x(1:nX) + x(2:nX+1))/2;  y_avg = (y(1:nY) + y(2:nY+1))/2;
            plot_diagonal_profile(x_avg, y_avg, reshape(v_vec, nX, nY), ...
                @(s) solveBurgers2(alpha, beta, 2*t_final, s), ...
                fullfile(pth, '2D', sprintf('rect_ls_burgers_t%s_N%d.pdf', strrep(num2str(t_final), '.', 'p'), nX)));
        end
    end
    T1 = order_table(DX,n_array, L1Error, "L1 Error");
    T2 = order_table(DX,n_array, L2_Error_quadratic, "L2 Error");
    T3 = order_table(DX,n_array, LinfError, "Linf Error");


    
%%%%%%Start of functions for multires%%%%%%%%%%%%

function u_hold = L(v_vec,W, quad_vals,quad_weights,len_vec,sten_array,coord_array,N,f,g,a)
    quad_res = weno_reconstruction(W,sten_array,coord_array,quad_vals,v_vec,N); 
    flux_array = find_fluxes(quad_res,sten_array,N,f,g,a);

    u_hold = NaN(N,1);
    for i = 1:N
        x_m = coord_array(i,1); x_p = coord_array(i,2); 
        y_m = coord_array(i,3); y_p = coord_array(i,4);
        sum = 0;
        for e = 1:4
            weights = quad_weights(i,e,:);
            weights = [0.5,0.5];
            len = len_vec(i,e);
            sum = sum + len*(weights(1)*flux_array(i,e,1) + weights(2)*flux_array(i,e,2));
        end
        area = (x_p - x_m)*(y_p - y_m);
        u_hold(i) = -(1/area)*sum;
    end
end

function quad_res = weno_reconstruction(W,sten_array,coord_array,quad_vals,v_vec,N)
   quad_res = NaN(N,4,2);
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

        u_cell = @(x,y) C(1) + C(2)*x + C(3)*y+ C(4)*x.*y + C(5)*x.^2 + C(6)*y.^2;

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
    a = 1;
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


  function [quad_vals, quad_weights, len_vec] = get_quad_info(coord_array,N)

  quad_vals = NaN(N,4,2); quad_weights = NaN(N,4,2);
        len_vec = NaN(N,4);
        for i = 1:N
          x_m = coord_array(i,1);  x_p = coord_array(i,2); 
          y_m = coord_array(i,3);  y_p = coord_array(i,4);
          
          dy = y_p - y_m;
          dx = x_p - x_m;

          %Edge 1 & 3
          [vals,weights] = lgwt(2,y_m,y_p); 
          quad_vals(i,1,:) = vals;
          quad_weights(i,1,:) = weights;

          quad_vals(i,3,:) = vals;
          quad_weights(i,3,:) = weights;

          len_vec(i,1) = dy;
          len_vec(i,3) = dy;

          %Edge 2 & 4
          [vals,weights] = lgwt(2,x_m,x_p); 
          quad_vals(i,2,:) = vals;
          quad_weights(i,2,:) = weights;

          quad_vals(i,4,:) = vals;
          quad_weights(i,4,:) = weights;

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
