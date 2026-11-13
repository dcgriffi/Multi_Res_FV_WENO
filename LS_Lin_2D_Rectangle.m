%LS_Lin_2D_Rectangle.m
%   Constrained least squares reconstruction of a linear polynomial
%   (k = 3 coefficients: 1, x, y) from cell averages on a rectangular grid.
%   Proof of concept for one one-sided direction of a linear stage: the
%   stencil is the target cell and its east and north neighbors (m = 3), so
%   the polynomial interpolates the three cell averages. A full linear stage
%   would repeat this for each one-sided direction.
%
%   test = "example"     : non-uniform grid on [0,1]^2. For one target cell,
%                          prints the difference of the cell averages of p
%                          and f over the whole grid and saves a figure of it
%                          to Results/LS_2D/ls_rect_lin_example.pdf.
%   test = "convergence" : L2 error of the reconstruction of sin(2*pi*(x+y))
%                          on uniform periodic grids of [0,1]^2. Expected
%                          order: 2.
%
%   The least squares system is built in local coordinates centered at the
%   target cell and scaled by sqrt(area), then mapped back to 1, x, y.
%   Requires lgwt.m.

clear; close all;

test = "convergence";   % "example" or "convergence"

k = 3; m = 3;            % linear polynomial; target, east and north cells
u_initial = @(x,y) sin(2*pi*(x+y));     % convergence test
f_example = @(x,y) sin(2*pi*(x+y));     % worked example

if test == "example"
    %Non-uniform grid on [0,1]^2
    x = [0; 0.5; 0.9; 1.4; 2.3; 3.2; 4]/4;
    y = [0, 0.6, 0.7, 1.2, 2.7, 3.5, 4.1, 5.2]/5.2;
    dx = x(2:end) - x(1:end-1); dy = y(2:end) - y(1:end-1);
    nX = length(x)-1; nY = length(y)-1; N = nX*nY;
    [ind_array,sten_array, coord_array, sten_coord_array] = get_rectangle_mesh_info(nX,nY,x,y,dx,dy,N);
    v_vec = build_cell_avg_vec(coord_array,f_example,N);
    i = ind_array(3,4);                          % interior target cell
    linear_stencil = [1,2,4];                    % target, east, north
    stencil = sten_array(i,:);
    W = rectangle_linear_condensed_constrained_weight_matrix(k,m,i,sten_coord_array,coord_array);
    C = W*[v_vec(stencil(linear_stencil)); -v_vec(i)];
    p = @(xx,yy) C(1) + C(2)*xx + C(3)*yy;
    used_cells = stencil(linear_stencil);
    pbar = zeros(N,1);
    for j = 1:N
        c = coord_array(j,:);
        pbar(j) = Gauss_quad_2D(c(1),c(2),c(3),c(4),p)/((c(2)-c(1))*(c(4)-c(3)));
    end
    diff_array = reshape(pbar - v_vec, [nX,nY]);    % (x index, y index)
    disp('Difference of cell averages (p - f), rows = x index, columns = y index:')
    disp(diff_array)
    in_stencil = false(nX,nY);
    for j = used_cells, in_stencil(ind_array == j) = true; end
    pth = fullfile(pwd, 'Results', 'LS_2D');
    if ~exist(pth, 'dir'), mkdir(pth); end
    plot_stencil_grid(diff_array, ind_array, i, in_stencil, fullfile(pth, 'ls_rect_lin_example.pdf'));
    return
end

n_array = [5,10,20,40,80,160];
DX = zeros(length(n_array),1); L2_Error_linear = zeros(length(n_array),1);
xL = 0; xR = 1;
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
        v_array = reshape(v_vec, [nX,nY]);

       

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


        global_L2_linear = global_L2_error(C,u_initial,N,coord_array,k);
        L2_Error_linear(kk)= (global_L2_linear)^(1/2);

    end

    T2 = order_table(DX,n_array, L2_Error_linear, "L2 Error");







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

function result = Gauss_quad_2D(xlow,xhigh,ylow,yhigh,interior)
        [x_vals,x_weights] = lgwt(5,xlow, xhigh);
        [y_vals, y_weights] = lgwt(5,ylow, yhigh);

        integrand = interior(x_vals,y_vals');
        result =  sum(sum(integrand.*x_weights.*y_weights'));
end


 function W = rectangle_linear_condensed_constrained_weight_matrix(k,m,i,sten_coord_array,coord_array)
        %Center the least squares on cell i and scale by sqrt(area). Building
        %the system in absolute coordinates makes A'A badly conditioned.
        xc = (coord_array(i,1) + coord_array(i,2))/2; yc = (coord_array(i,3) + coord_array(i,4))/2;
        s = sqrt((coord_array(i,2) - coord_array(i,1))*(coord_array(i,4) - coord_array(i,3)));


        %Build matrix for least squares
        A = zeros(m,k);
        linear_stencil = [1,2,4,6,8];
        for j = 1:m
            cell = linear_stencil(j);
            x_m = (sten_coord_array(i,cell,1)-xc)/s; x_p = (sten_coord_array(i,cell,2)-xc)/s; y_m = (sten_coord_array(i,cell,3)-yc)/s; y_p = (sten_coord_array(i,cell,4)-yc)/s;
            area = (x_p - x_m)*(y_p - y_m);
            A(j,1) = 1/area*(x_p - x_m)*(y_p - y_m);
            A(j,2) = 1/area*(x_p^2 - x_m^2)*(y_p-y_m)/2;
            A(j,3) = 1/area*(y_p^2 - y_m^2)*(x_p-x_m)/2;
        end
        
        %Build L matrix
        L = zeros(1,k);
        cell = i;
        x_m = (coord_array(cell,1)-xc)/s; x_p = (coord_array(cell,2)-xc)/s; y_m = (coord_array(cell,3)-yc)/s; y_p = (coord_array(cell,4)-yc)/s;
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

        %Map coefficients of 1, z, w (z = (x-xc)/s, w = (y-yc)/s) back to
        %coefficients of 1, x, y, so the rest of the code is unchanged.
        W = [1, -xc/s, -yc/s; 0, 1/s, 0; 0, 0, 1/s]*W;


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

function plot_stencil_grid(diff_array, ind_array, i, in_stencil, fname)
    %Grid of cells in index layout (x to the right, y up) with the value of
    %diff_array in each cell. Stencil cells are shaded light gray and
    %outlined in red, and the target cell i is shaded darker. No title;
    %details go in the LaTeX caption.
    [nX, nY] = size(diff_array);
    [ti, tj] = find(ind_array == i);
    fig = figure('Units','inches','Position',[1 1 8.0 5.6],'Color','w'); hold on
    for ix = 1:nX
        for jy = 1:nY
            if ix == ti && jy == tj
                col = [0.75 0.75 0.75];
            elseif in_stencil(ix,jy)
                col = [0.92 0.92 0.92];
            else
                col = [1 1 1];
            end
            rectangle('Position',[ix-1, jy-1, 1, 1], 'FaceColor', col, 'EdgeColor', [0.5 0.5 0.5]);
            text(ix-0.5, jy-0.5, sprintf('%.1e', diff_array(ix,jy)), ...
                'HorizontalAlignment','center', 'FontSize', 8);
        end
    end
    for ix = 1:nX
        for jy = 1:nY
            if in_stencil(ix,jy)
                %Outline the outer edges of the stencil in red
                if ix == 1  || ~in_stencil(ix-1,jy), plot([ix-1 ix-1],[jy-1 jy],'r-','LineWidth',2); end
                if ix == nX || ~in_stencil(ix+1,jy), plot([ix ix],[jy-1 jy],'r-','LineWidth',2); end
                if jy == 1  || ~in_stencil(ix,jy-1), plot([ix-1 ix],[jy-1 jy-1],'r-','LineWidth',2); end
                if jy == nY || ~in_stencil(ix,jy+1), plot([ix-1 ix],[jy jy],'r-','LineWidth',2); end
            end
        end
    end
    axis equal; xlim([0 nX]); ylim([0 nY]);
    xticks(0.5:1:nX-0.5); xticklabels(arrayfun(@(q) sprintf('%d', q), 1:nX, 'UniformOutput', false));
    yticks(0.5:1:nY-0.5); yticklabels(arrayfun(@(q) sprintf('%d', q), 1:nY, 'UniformOutput', false));
    xlabel('Cell index in x'); ylabel('Cell index in y');
    set(gca, 'FontSize', 10, 'TickLength', [0 0], 'Box', 'on');
    exportgraphics(fig, fname, 'ContentType', 'vector');
end
