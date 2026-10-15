%Convergence test for LS linear polynomials with uniform rectangular mesh (constrained) 


clear all; close all;

pth = '';


%Test for second order approximation

n_array=[5,10,20,40,80,160]; 

%Preallocate arrays for memory
DX = zeros(length(n_array),1); dt = zeros(length(n_array),1); 
L2_Error_quadratic = zeros(length(n_array),1); L1Error = zeros(length(n_array),1); 
L2_Error_linear = zeros(length(n_array),1);
LinfError= zeros(length(n_array),1); err = zeros(length(n_array),1); 

%In each cell, obtain polynomial approximation and store coefficients
%as we did for DG

xL = 0; xR = 1;
yL = xL; yR = xR;

%Cycle through all ranges of dx
   p = 3;
   u_initial= @(x,y) sin(2*pi*(x+y));
   %u_initial= @(x,y) 3 + 4*x + 3*y -4*x.*y + 7*x.^2 + 10*y.^2;
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
        v_array = reshape(v_vec, [nX,nY]);

       
        %Second order approximation
        k = 3; m = 3; 

        %Obtain weight matrix for coefficents
        W = zeros(k, m+1, N);
        for i = 1:N
            W(:,:,i) = rectangle_linear_condensed_constrained_weight_matrix(k,m,i,sten_coord_array,coord_array);
        end
    
        C = zeros(N,k); linear_stencil = [1,2,4,6,8];
        %Build coefficient matrix using weight matrix
        for i = 1:N
            stencil = sten_array(i,:);
            d = zeros(m+1,1);
            for j = 1:m
                index = linear_stencil(j);
                cell = stencil(index);
                d(j) = v_vec(cell);
            end
            d(m+1) = -v_vec(i);
            C(i,:) = W(:,:,i)*d;
        end


        %Prep for visualization
        x_dense = (xL:h/10:xR)';
        y_dense = (yL:h/10:yR);
        u_ref_x = u_initial(x_dense, 0.001);
        u_ref_y = u_initial(0.001, y_dense);
        fig_num = kk;
        fig_num = fig_num+1; variable_name = "sin(x+y)"; 
        Method = "Least Squares Poly Approx"; 
        
        if nX == 10 || nX == 160
             plot_solution(variable_name,C,u_ref_x, u_ref_y, fig_num,Method,N,nX,nY,x,y,x_dense,y_dense,h,k)
            % saveas(fig_num, sprintf('%s/Least_Square_Rectangle/least_square_plot_x_linear_%d.jpg',pth,nX),'jpeg')
            % saveas(fig_num+1, sprintf('%s/Least_Square_Rectangle/least_square_plot_y_linear_%d.jpg',pth,nX),'jpeg')
        end
         
        global_L2_linear = global_L2_error(C,u_initial,N,coord_array,k);
        L2_Error_linear(kk)= (global_L2_linear)^(1/2);

        %Plot 2D solution for nX = 10 (quadratic, three points)\
        fig_num = fig_num + 1;
        x_dense = NaN(30,1);
        y_dense = NaN(1,30);
        u_approx_array = NaN(30,30);
        u_exact_array = NaN(30,30);

        if nX == 10
           plot_2D_rectangle_exact_approx(x,y,nX,nY,C,ind_array,u_initial,fig_num,h,pth)
        end


        fig_num = fig_num + 1;     
                 
    end

    p = 2;
    fig_num = kk + 20; err_name = "L2 Error";
    T2 = order_table(DX,n_array, L2_Error_linear, err_name);
    order_plot(L2_Error_linear,fig_num,Method,p,DX,err_name, k);

   % saveas(fig_num, sprintf('%s/Least_Square_Rectangle/constrained_least_square_lin_convergence_plot.jpg',pth),'jpeg')







function order_plot(L2Error,fig_num,Method,p,dx,err_name, k)
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


function global_L2_Error = global_L2_error(C,exact_sol,N,coord_array,k)
    global_L2_Error = 0;
    for j = 1:N
        cell = j;
        %Obtain approx function at each cell (change if increase k)
        if k == 6
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
        elseif k == 3
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y;
        end
        x_m = coord_array(j,1); x_p = coord_array(j,2); y_m = coord_array(j,3); y_p = coord_array(j,4);

        %Obtain local error estimate with quadrature
        interior= @(x,y) (u_cell(x,y) - exact_sol(x,y)).^2;
        local_error = Gauss_quad_2D(x_m,x_p,y_m,y_p,interior); %Obtain integral via  Legendre-Gauss quadrature
        global_L2_Error = global_L2_Error + local_error;
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
        if k == 6
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
        elseif k == 3
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y;
        end
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
        
        if k == 6
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y+ C(cell ,4)*x.*y + C(cell ,5)*x.^2 + C(cell ,6)*y.^2;
        elseif k == 3
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y;
        end

        %Plot 5 points in the interval  + 
        int_points = (y(i):h/5:y(i+1)); %This enumerates cells in x direction, midslice tells which cells in y. 
        plot(int_points,u_cell(0.001, int_points))
    end
    title(sprintf("Slice of %s at x = 0 from %s with nX = ", [variable_name Method]) + nY)
    legend('Exact','Approximate')
end

function result = Gauss_quad_2D(xlow,xhigh,ylow,yhigh,interior)
        [x_vals,x_weights] = lgwt(5,xlow, xhigh);
        [y_vals, y_weights] = lgwt(5,ylow, yhigh);

        integrand = interior(x_vals,y_vals');
        result =  sum(sum(integrand.*x_weights.*y_weights'));
end


 function W = rectangle_linear_condensed_constrained_weight_matrix(k,m,i,sten_coord_array,coord_array)

        %Build matrix for least squares
        A = zeros(m,k);
        linear_stencil = [1,2,4,6,8];
        for j = 1:m
            cell = linear_stencil(j);
            x_m = sten_coord_array(i,cell,1); x_p = sten_coord_array(i,cell,2); y_m = sten_coord_array(i,cell,3); y_p = sten_coord_array(i,cell,4);
            area = (x_p - x_m)*(y_p - y_m);
            A(j,1) = 1/area*(x_p - x_m)*(y_p - y_m);
            A(j,2) = 1/area*(x_p^2 - x_m^2)*(y_p-y_m)/2;
            A(j,3) = 1/area*(y_p^2 - y_m^2)*(x_p-x_m)/2;
        end
        
        %Build L matrix
        L = zeros(1,k);
        cell = i;
        x_m = coord_array(cell,1); x_p = coord_array(cell,2); y_m = coord_array(cell,3); y_p = coord_array(cell,4);
        area = (x_p - x_m)*(y_p - y_m);
        L(1) = 1/area*(x_p - x_m)*(y_p - y_m);
        L(2) = 1/area*(x_p^2 - x_m^2)*(y_p-y_m)/2;
        L(3) = 1/area*(y_p^2 - y_m^2)*(x_p-x_m)/2;
        
        
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
                 u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y;
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
    title(sprintf("2D Plot of constrained least square for sin(x + y)"))
    %saveas(fig_num, sprintf('%s/Least_Square_Rectangle/least_square_plot_2D_linear_%d.jpg',pth,nX),'jpeg')
    
    fig_num = fig_num + 1;
    figure(fig_num)
    surf(X,Y,u_exact_array);
    title(sprintf("2D Plot of exact solution for sin(x + y)"))
    %saveas(fig_num, sprintf('%s/Least_Square_Rectangle/exact_plot_2D_linear_%d.jpg',pth,nX),'jpeg')
 end


 function difference = test_cell_avg_difference_rectangle(coord_array,nX,nY,N,C,v_array);
 %Check the difference between the function and cell averages for
        %each cell (it works)
        
        cell_avg = zeros(N,1);
        for i = 1:N
            cell = i;
            %Compute cell average using polynomial c
            u_cell = @(x,y) C(cell ,1) + C(cell ,2)*x + C(cell ,3)*y;
            x_m = coord_array(i,1); x_p = coord_array(i,2); y_m = coord_array(i,3); y_p = coord_array(i,4);
            area = (x_p - x_m)*(y_p- y_m);
            %Obtain local error estimate with quadrature
            interior= @(x,y) u_cell(x,y);
            cell_avg(i) = Gauss_quad_2D(x_m,x_p,y_m,y_p,interior)/area;

        end
        

        cell_avg_array = reshape(cell_avg, [nX,nY]);

        difference = abs(v_array - cell_avg_array);
 end

            