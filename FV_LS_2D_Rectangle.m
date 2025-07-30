%Convergence test for FV LS 3rd order polynomials with uniform rectangular mesh (constrained) 
    %%% Run pde = linear for order or pde = burgers for behavior around shock

clear all; close all;

%pth = '/Users/chloegriffin/Brown Dropbox/Chloe Griffin/Research 2024/AWI Codes/Results';
pth = '';
%Uncomment the desired time
%t_final = 2/2; %After shock
t_final = 0.3/2;

avoid_shock = true;

p = 3;

%pde = "burgers";
pde = "linear";
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

n_array=[5,10,20,40,80]; 

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

        x_dense = (xL:h/10:xR)';
        y_dense = (yL:h/10:yR);
        u_ref_x = u_ref(x_dense, 0.001);
        u_ref_y = u_ref(0.001, y_dense);
        fig_num = kk;
        fig_num = fig_num+1; variable_name = "sin(x+y)"; 
        Method = "Least Squares Poly Approx"; 

        %Updating this to find error using cell average (Old way)
        C = zeros(N,k);
        %Build coefficient matrix using weight matrix
        for i = 1:N
            stencil = sten_array(i,:);
            d = zeros(m+1,1);
            for j = 1:m
                cell = stencil(j);
                d(j) = v_vec(cell);
            end
            d(m+1) = -v_vec(i);
        
            C(i,:) = W(:,:,i)*d;
        end
        
       

        if nX == 80
             plot_solution(variable_name,C,u_ref_x, u_ref_y, fig_num,Method,N,nX,nY,x,y,x_dense,y_dense,h,k)
            % saveas(fig_num, sprintf('%s/FV_Least_Square_Rectangle/least_square_plot_x_quadratic_%d_%d.jpg',pth,nX,t_final),'jpeg')
            % saveas(fig_num+1, sprintf('%s/FV_Least_Square_Rectangle/least_square_plot_y_quadratic_%d_%d.jpg',pth,nX,t_final),'jpeg')
             fig_num = fig_num + 3;
        end
    
        u_exact_avg = build_cell_avg_vec(coord_array,u_ref,N);
        
        if avoid_shock == false
        %Compute uRef
        uDif = v_vec - u_exact_avg;
        L2_Error_quadratic(kk) = ((1/N)*sum(abs(uDif).^2))^(1/2);
        end

        if avoid_shock == true
        [cell_count, global_L2_quadratic] = global_L2_error(N,coord_array,t_final,alpha,u_exact_avg,v_vec);
        L2_Error_quadratic(kk)= ((1/cell_count)*(global_L2_quadratic))^(1/2);
        end

        %Plot 2D solution for nX = 10 (quadratic, three points)
        x_dense = NaN(30,1);
        y_dense = NaN(1,30);
        u_approx_array = NaN(30,30);
        u_exact_array = NaN(30,30);

        if nX == 40
           plot_2D_rectangle_exact_approx(x,y,nX,nY,C,ind_array,u_ref,fig_num,h,pth)
        end

        fig_num = fig_num + 1;
    end
    format shorte
    p = 3;
    fig_num = fig_num+2; err_name = "L2 Error";
    T1 = order_table(DX,n_array, L2_Error_quadratic, err_name);
    order_plot(L2_Error_quadratic,fig_num,Method,p,DX,err_name, k);
    %saveas(fig_num, sprintf('%s/FV_Least_Square_Rectangle/ls_convergence_plot_%d.jpg',pth, t_final),'jpeg')


    
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
function order_plot(L2Error,fig_num,Method,p,dx,err_name, k)
    figure(fig_num)
    hold on;
    plot(log10(dx),log10(L2Error),'.-','LineWidth',2,'MarkerSize',20)
    plot(log10(dx),p*log10(dx),'k','LineWidth',3) %Reference Slope = 3
    xlabel('Log_{10}(\Deltax)')
    ylabel('Log_{10}(Error)')
    title(sprintf("%s Convergence Plot using %s with K = %s", [err_name Method k]))
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


% function global_L2_Error = global_L2_error(C,exact_sol,N,coord_array,avoid_shock,t_final, alpha)
function [cell_count, global_L2_Error] = global_L2_error(N,coord_array,t_final,alpha,u_exact_avg,v_vec)
    global_L2_Error = 0;
  
    xShock = 2*pi + alpha * t_final;
    diag_shock1 = xShock - pi;
    diag_shock2 = xShock + pi;

    low_1 = diag_shock1 - 1;
    low_2 = diag_shock2 - 1;

    high_1 = diag_shock1 + 1;
    high_2 = diag_shock2 + 1;
     
    cell_count = 0;
    for j = 1:N
        x_m = coord_array(j,1); x_p = coord_array(j,2); y_m = coord_array(j,3); y_p = coord_array(j,4);
        
        if ((x_p + y_p) <= low_1) || (((x_m + y_m) >= high_1) && ((x_p + y_p) <= low_2)) || ((x_m + y_m) >= high_2)
            %Obtain local error estimate with quadrature
            local_error = abs(v_vec(j) - u_exact_avg(j))^2; %Obtain integral via  Legendre-Gauss quadrature
            global_L2_Error = global_L2_Error + local_error;

            cell_count = cell_count + 1;
        end
    end
end

function plot_solution(variable_name,C,u_ref_x, u_ref_y, fig_num,Method,N,nX,nY,x,y,x_dense,y_dense,h,k)

%Update to do one plot in x, fixing y and one plot in y, fixing x = 0
    %Plot Solution (fix y = 0)
    figure(fig_num)
    plot(x_dense,u_ref_x) %Plot exact (need u_ref inside)
    hold on

    for i = 1:nX
        cell = i;
        u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
       
        %Plot 5 points in the interval  + 
        int_points = (x(i):h/5:x(i+1));%This enumerates cells in x direction, midslice tells which cells in y. 
        plot(int_points,u_cell(int_points, 0.001))
    end
    title(sprintf("Slice of %s at y = 0 from %s with nX = ", [variable_name Method]) + nX)
    legend('Exact','Approximate')

    %Plot Solution (fix x = 0)
    figure(fig_num+1)
    plot(y_dense,u_ref_y) %Plot exact (need u_ref inside)
    hold on
    for i = 1:nY
        cell = 1 + (i - 1)*nX;
        
        u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
       

        %Plot 5 points in the interval  + 
        int_points = (y(i):h/5:y(i+1)); %This enumerates cells in x direction, midslice tells which cells in y. 
        plot(int_points,u_cell(0.001, int_points))
    end
    title(sprintf("Slice of %s at x = 0 from %s with nX = ", [variable_name Method]) + nY)
    legend('Exact','Approximate')
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

         %Build matrix for least squares
        A = zeros(m,k);
        for j = 1:m
            x_m = sten_coord_array(i,j,1); x_p = sten_coord_array(i,j,2); y_m = sten_coord_array(i,j,3); y_p = sten_coord_array(i,j,4);
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
        x_m = coord_array(cell,1); x_p = coord_array(cell,2); y_m = coord_array(cell,3); y_p = coord_array(cell,4);
        area = (x_p - x_m)*(y_p - y_m);
        L(1) = 1/area*(x_p - x_m)*(y_p - y_m);
        L(2) = 1/area*(x_p^2 - x_m^2)*(y_p-y_m)/2;
        L(3) = 1/area*(y_p^2 - y_m^2)*(x_p-x_m)/2;
        L(4) = 1/area*(x_p^2 - x_m^2)*(y_p^2 - y_m^2)/4;
        L(5) = 1/area*(x_p^3 - x_m^3)*(y_p-y_m)/3;
        L(6) = 1/area*(y_p^3 - y_m^3)*(x_p-x_m)/3;
        
        
        %Try again with QR Factorization
        
        [Q,R] = qr(A);
        W = R\(Q');    %(A'*A)\A'
        W = (A'*A)\(A');
        %P = R\(Q'*L'); %(A'*A)\(L');
        P =(A'*A)\(L');
        
        mu = -L*P;
        E = W + (P*L*W)/mu;
        F = P/mu;
        
        W = [E, F];


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


 function plot_2D_rectangle_exact_approx(x,y,nX,nY,C,ind_array,u_initial,fig_num,h,pth)
    fig_num = fig_num + 1;
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
    title(sprintf("2D Plot of constrained least square for 1/3 + 2/3sin((x + y)/2"))
   % saveas(fig_num, sprintf('%s/FV_Least_Square_Rectangle//least_square_plot_2D_quadratic_%d.jpg',pth,nX),'jpeg')
    
    fig_num = fig_num + 1;
    figure(fig_num)
    surf(X,Y,u_exact_array);
    title(sprintf("2D Plot of exact solution for 1/3 + 2/3sin((x + y)/2)"))
    %saveas(fig_num, sprintf('%s/FV_Least_Square_Rectangle/exact_plot_2D_quadratic_%d.jpg',pth,nX),'jpeg')
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

 