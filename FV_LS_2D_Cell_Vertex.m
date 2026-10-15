%Convergence test for FV LS 3rd order polynomials with uniform cell-vertex mesh (constrained) 
    %%% Run pde = linear for order or pde = burgers for behavior around shock

clear all; close all;


    %Stored Info Guide: 
    % n2d: Number of nodes/control volumes/ triangle vertices (nx x ny)
    % e2d: Number of triangles

    % xcoord: (1xn2d) x-coordinates of control volume center/triangle vertices
    % ycoord  (1xn2d) y-coordinates of control volume center/triangle vertices
    
    % tri: (3xe2d) Tracks three vertices associated to each triangle (used
            % for plotting mesh)
    % nod_in_e2d: (6xn2d) Tracks triangles for each node (unordered)
    % nghbr_n2d:  (7xn2d) Tracks neighboring nodes for each node (unordered
            % for now)
    % ng_num: (1xn2d) Tracks num of neighbors for each node 
    % n_num:  (1xn2d) Tracks num of triagles for each node
    
    % edge_num/ed2d: Number of edges (not sure why both)
    % edges: (2xedge_num) Tracks vertices for each edge
    % edge_tri : (2xedge_num) Tracks triangle/triangles associated to each edge
    % elem_edges: (3xe2d) Tracks edges associated to each triangle (numbered as in
    %      edges)

    % edge_dxdy : (2xedge_num) Gives dx and dy for each edge
    % edge_cross_dxdy: (4xedge_num) Gives cross-distances for the edges

    %enodes : tri(:,n) Nodes of triangle for a particular control volume

%Uncomment the desired time
%t_final = 2/2; %After shock
t_final = 0.3/2; %Before shock
%t_final = 1e-10;


p = 3;

pde = "burgers";
%pde = "linear";
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


n_array = [9,17,33,65,129,257];
%n_array = [10,20,40,80,160,320];
n_array = [20,40,80,160];
%n_array = [80];




Method = "Least Squares Poly Approx"; p = 3;
    

%Preallocate arrays for memory
DX = zeros(length(n_array),1); dt = zeros(length(n_array),1); 
L2_Error_quadratic = zeros(length(n_array),1); L1Error = zeros(length(n_array),1); 
L2_Error_linear = zeros(length(n_array),1);
LinfError= zeros(length(n_array),1); err = zeros(length(n_array),1); 

n_set_array = n_array;

xL = 0; xR = 2*pi;
yL = 0; yR = 2*pi;

   %Cycle through all ranges of dx
   for kk=1:length(n_array)

        %Build spatial array
         n = n_set_array(kk); 
         h=(xR-xL)/n;
         x = (xL:h:xR);
         y = (yL:h:yR); 
        
         dx  = x(2:end) - x(1:end-1);
         dy  = y(2:end) - y(1:end-1);

         %Set up ghost points

         x = [x(1) - dx(end) - dx(end - 1) - dx(end - 2),x(1) - dx(end) - dx(end - 1) , x(1) - dx(end), x, x(end) + h, x(end) + dx(1) + dx(2), x(end) + dx(1) + dx(2)+dx(3)];
         y = [y(1) - dy(end) - dy(end - 1) - dy(end - 2),y(1) - dy(end) - dy(end - 1),  y(1) - dx(end), y, y(end) + h, y(end) + dy(1) + dy(2), y(end) + dy(1) + dy(2)+dy(3)];
       
         N = n + 7;

         DX(kk)=(xR-xL)/(n);
       
         [tri, nodnum, xcoord, ycoord, e2d, n2d] = gen_reg_tri_mesh(x,y);

         [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d, edge_tri,edge_cross_dxdy] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d);
        
         bound_tri_vec = 1:N*N;
        
         bound_tri_array = reshape(bound_tri_vec, N,N);
         inner_v_array = bound_tri_array(2:end - 1, 2:end - 1);%Array for all cell averages

         inner_v_array_2 = inner_v_array(2:end-1, 2:end-1); %Array for all ghost reconstruction
         ind_array = inner_v_array(3:end - 2, 3:end - 2);
        
         inner_v_vec = reshape(inner_v_array,1,(N-2)*(N-2));
         inner_v_vec_2 = reshape(inner_v_array_2, 1,(N-4)*(N-4));

         ind_vec = reshape(ind_array, 1,(n+1)*(n+1)); %or (N-6)*(N-6)
         n_steps = (n+1)*(n+1);
         inner_steps = (N-4)*(N-4);
        
         cell_avg = build_cell_avg_vec(inner_v_vec, u_initial, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
         %boundary_check = apply_bc(cell_avg,N);

         cell_avg_array = reshape(cell_avg, N, N);
      
         %Third order approximation
         k = 6; m = 7;
 

          W = NaN(k, m+1, N*N);
          for ii = 1:inner_steps
            i = inner_v_vec_2(ii);
            W(:,:,i) = hexagon_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area);
          end

         u = cell_avg; %(Note: cell_avg = NxN - Only update (N-4xN-4),(N-2xN-2) with periodic boundary

        %Change both functions into plottable array
        t = 0;
        while t < t_final     
            %Choose appropriate time step
             h = max(dx);
             a = max(abs(u));
             dt = (h)/(4*a); 
             %dt = h*0.07;

            if t_final - t < dt
                 dt = t_final - t;
            end
    
            %RK3
            u_1 = u + dt*L(u,n_steps, N, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,edge_cross_dxdy);
            u_1 = apply_bc(u_1,N);

            u_2 = 3/4*u +1/4*u_1 + 1/4*dt*L(u_1,n_steps, N, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,edge_cross_dxdy);
            u_2 = apply_bc(u_2,N);

            u_next = 1/3*u + 2/3*u_2 + 2/3*dt*L(u_2,n_steps, N, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,edge_cross_dxdy);
            u_next = apply_bc(u_next,N);
    
            u = u_next;
            t = t + dt;
        end
        
        %Prep for visualization

        if pde == "burgers"
            %u_ref= @(x,y) solveBurgers2(alpha, beta, t_final,(x+y)./2);
            u_ref= @(x,y) solveBurgers2(alpha, beta, 2*t_final,(x+y));
             %u_ref= @(x,y) solveBurgers2(alpha, beta, 2*t_final,(x+y));

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
       
        %Visualization: Need NxN with NaNs for visualization, cut off for
        %error

        if n == n_set_array(end)
             %Plot exact cell averages
             exact_avg_vec = build_cell_avg_vec(inner_v_vec, u_ref, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
             fig_num = 3;
             exact_avg_array = reshape(exact_avg_vec,N,N);
             figure(fig_num)
             trimesh(tri', xcoord, ycoord, exact_avg_array)
             title("Plot of Exact Cell Avg at t=" + t_final + " with nX = " + N)
             %saveas(fig_num, sprintf('%s/FV_Least_Square_Hexagon/least_square_plot_exact_%d_%d.jpg',pth,N,t_final),'jpeg')



             %Plot approx cell averages
             fig_num = 4;
             approx_avg_array = reshape(u,N,N);
             figure(fig_num)
             trimesh(tri', xcoord, ycoord, approx_avg_array)
             title("Plot of Approx Cell Avg at t=" + t_final + " with nX = " + N)
             %saveas(fig_num, sprintf('%s/FV_Least_Square_Hexagon/least_square_plot_exact_%d_%d.jpg',pth,N,t_final),'jpeg')

             fig_num = fig_num + 3;
        end
    
        exact_avg_vec = build_cell_avg_vec(inner_v_vec, u_ref, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
        
        exact_avg_array = reshape(exact_avg_vec,N,N);
        approx_avg_array = reshape(u,N,N);

       
        %Compute Error
        uDif = approx_avg_array(4:end - 3, 4:end -3) - exact_avg_array(4:end - 3, 4:end - 3);

        L2_Error_quadratic(kk) = ((1/(n_steps))*sum(sum(abs(uDif).^2)))^(1/2);
   
   end
   
    format shorte
    p = 3;
    err_name = "L2 Error";
    T1 = order_table(DX,n_array, L2_Error_quadratic, err_name);
    order_plot(L2_Error_quadratic, Method,p,DX,err_name, k);
    %saveas(6, sprintf('%s/FV_Least_Square_Hexagon/ls_convergence_plot_%d.jpg',pth, t_final),'jpeg')

    % 
    % %Test periodic boundary on exact solution

    % u_1 = apply_bc(exact_avg_vec,N);
    % u_1 = reshape(u_1,N,N);
    % shock_check = u_1 - exact_avg_array
    % 
    % x_array = xcoord(bound_tri_array);
    % x_vec = reshape(x_array,1,N*N);
    % 
    % u_1 = apply_bc(x_vec,N);
    % u_1 = reshape(u_1,N,N);
    % x_check = u_1 - x_array
    % 
    % 
    % y_array = ycoord(bound_tri_array);
    % y_vec = reshape(y_array,1,N*N);
    % 
    % u_1 = apply_bc(y_vec,N);
    % u_1 = reshape(u_1,N,N);
    % y_check = u_1 - y_array
    % 
    % area_array = reshape(area,N,N);
    % area_vec = reshape(area_array,1,N*N);
    % 
    % u_1 = apply_bc(area_vec,N);
    % u_1 = reshape(u_1,N,N);
    % area_check = u_1 - area_array





%%%%%%%%%%%Start of new functions%%%%%%%%%%%%%%%%%%%%%%%%%%


function u_1 = apply_bc(u_1,N)
u_1 = reshape(u_1,N,N);

%Pad twice periodically
u_1(3,:) = u_1(end - 4,:);
u_1(end-2,:) = u_1(5,:);
u_1(:,3) = u_1(:,end - 4);
u_1(:,end-2) = u_1(:,5);

u_1(2,:) = u_1(end - 5,:);
u_1(end-1,:) = u_1(6,:);
u_1(:,2) = u_1(:,end - 5);
u_1(:,end-1) = u_1(:,6);

u_1 = reshape(u_1,1,N*N);
end

function u_hold = L(u,n_steps, N, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,edge_cross_dxdy)

    inner_sum = find_inner_sum(u,W,ed2d,xcoord,ycoord,edges,edge_tri,tri,f,g,N,nghbr_n2d,inner_v_vec_2,inner_steps,edge_cross_dxdy);

    u_hold = NaN(1,N*N);

    %This can easily be vectorized, keeping this way for Fortran use
    for ii = 1:n_steps
        i = ind_vec(ii);
        sum =0;
        for e = 1:12
            sum = sum + inner_sum(e,i);
        end
        area_i = area(i);
        u_hold(i) = -(1/area_i)*sum;
    end
end


function inner_sum = find_inner_sum(u,W,ed2d,xcoord,ycoord,edges,edge_tri,tri,f,g,N,nghbr_n2d,inner_v_vec_2,inner_steps,edge_cross_dxdy)

inner_sum = zeros(12,N*N);
counter = ones(1,N*N); %Tracks which index to fill flux on edge (remove later)

%For LS reconstruction, store coefficients for all i 
C = weno_reconstruction(u,W,nghbr_n2d,N,inner_v_vec_2,inner_steps);

%Third order reconstruction
u_cell = @(x,y,i) C(1,i) + C(2,i)*x + C(3,i)*y+ C(4,i)*x.*y + C(5,i)*x.^2 + C(6,i)*y.^2;

%First order reconstruction
%u_cell = @(x,y,i) u(i);


d = @(v1,v2) sqrt((v1(1) - v2(1))^2 + (v1(2) - v2(2))^2);
T = @(s,a,b) a*(1 - s)/2+ b*(1 + s)/2; %Transform from reference segment

%%%Next steps: Need to detect edges attached only to interior points

for e = 1:ed2d
    v1 = edges(1,e); v2 = edges(2,e);

    %Only want to compute results for interior points (where two triangles
    %adjacent) %Maybe can run until edge_in
    if (edge_tri(1,e) && edge_tri(2,e) ~= 0)
        cross_info = edge_cross_dxdy(:,e);
        quad_points = NaN(2,2);
        midpoint = [(xcoord(v1) + xcoord(v2))/2, (ycoord(v1) + ycoord(v2))/2];%Compute midpoint


        for side = 1:2
            %Compute centroid

            cross_vec = cross_info((2*side - 1): 2*side)';

            centroid = midpoint + cross_vec;

            %Get 2 quad points
            quad_points(:,1) = T(-1/sqrt(3),midpoint,centroid); quad_points(:,2) = T(1/sqrt(3),midpoint,centroid);
            
            %Get edge length
            edge_len = norm(cross_vec);
            %Get outward facing normal vector to side edge
            
            %Rotate 90 degrees/-90 degrees
            if side == 1
                 n_vec = [cross_vec(2); -cross_vec(1)]/edge_len;
            else 
                 n_vec = [-cross_vec(2); cross_vec(1)]/edge_len;
            end
            
            fv_1 = NaN(2,1); fv_2 = NaN(2,1);
            for q = 1:2
                %Use LS/multi-res for scalar reconstruction from either cell
                uq_1 = u_cell(quad_points(1,q),quad_points(2,q),v1);
                uq_2 = u_cell(quad_points(1,q),quad_points(2,q),v2);
            
                fv_1(q) = lf(uq_1,uq_2,n_vec,f,g); %u_m, u_p
                fv_2(q) = -fv_1(q);
            end
            inner_sum(counter(v1),v1) = edge_len*(1/2*(fv_1(1) + fv_1(2))); %Summing over quad points for edge
            inner_sum(counter(v2),v2) = edge_len*(1/2*(fv_2(1) + fv_2(2))); %weights = 1/2
            %inner_sum(counter(v1),v1) = edge_len*(1/2*(fv_1(1) + fv_1(2))); %Summing over quad points for edge
            %inner_sum(counter(v2),v2) = edge_len*(1/2*(fv_2(1) + fv_2(2))); %weights = 1/2
            counter(v1) = counter(v1) + 1;
            counter(v2) = counter(v2) + 1;
         end
     end
   
end


%Old approach without edge_cross_dxdy
% for e = 1:ed2d
%     v1 = edges(1,e); v2 = edges(2,e);
%     v1_coord = [xcoord(v1), ycoord(v1)]; v2_coord = [xcoord(v2), ycoord(v2)];
%     %Only want to compute results for interior points (where two triangles
%     %adjacent %Maybe can run until edge_in
%     if (edge_tri(1,e) && edge_tri(2,e) ~= 0)
%         quad_points = NaN(2,2);
%         midpoint = [(xcoord(v1) + xcoord(v2))/2, (ycoord(v1) + ycoord(v2))/2];%Compute midpoint
%         for side = 1:2
%             side_tri = edge_tri(side,e);
%             tri_vert = tri(:,side_tri);
%             %Compute centroid
%             centroid = [sum(xcoord(tri_vert))/3, sum(ycoord(tri_vert))/3];
%             %Get 2 quad points
%             quad_points(:,1) = T(-1/sqrt(3),midpoint,centroid); quad_points(:,2) = T(1/sqrt(3),midpoint,centroid);
%             %Get edge length
%             edge_len =d(centroid,midpoint); %sqrt((centroid(1) - midpoint(1))^2 + (centroid(2) - midpoint(2)^2));
%             %Get outward facing normal vector to side edge
%             n_vec = [-(centroid(2) - midpoint(2)),(centroid(1) - midpoint(1))]/(edge_len); %(-dy,dx)/edge_len
%             test = midpoint + n_vec;
%             fv_1 = NaN(2,1); fv_2 = NaN(2,1);
%             for q = 1:2
%                 %Use LS/multi-res for scalar reconstruction from either cell
%                 uq_1 = u_cell(quad_points(1,q),quad_points(2,q),v1);
%                 uq_2 = u_cell(quad_points(1,q),quad_points(2,q),v2);
%                 %Assign correct flux depending on direction of normal vector
%                 a = 1;
%                 if d(v1_coord,test) < d(v2_coord,test)
%                     fv_2(q) = lf(uq_2,uq_1,n_vec,f,g,a); %n_vec outward facing from 
%                     fv_1(q) = -fv_2(q);
%                 else 
%                     fv_1(q) = lf(uq_1,uq_2,n_vec,f,g,a); %u_m, u_p
%                     fv_2(q) = -fv_1(q);
%                 end
%             end
%             inner_sum(counter(v1),v1) = edge_len*(1/2*(fv_1(1) + fv_1(2))); %Summing over quad points for edge
%             inner_sum(counter(v2),v2) = edge_len*(1/2*(fv_2(1) + fv_2(2))); %weights = 1/2
%             counter(v1) = counter(v1) + 1;
%             counter(v2) = counter(v2) + 1;
%          end
%      end
% 
% end
end


function C  = weno_reconstruction(u,W,nghbr_n2d,N,inner_v_vec_2,inner_steps)
      k = 6; m = 7;
      C = NaN(k,N*N);
      %Build coefficient matrix using weight matrix
      for ii = 1:inner_steps
        i = inner_v_vec_2(ii);
        stencil = nghbr_n2d(:,i);
        d = zeros(m+1,1);
        for j = 1:m
            cell = stencil(j);
            d(j) = u(cell);
        end
        d(m+1) = -u(i);
    
        C(:,i) = W(:,:,i)*d;
      end
end

function fhat = lf(u_m,u_p,n_vec,f,g)
    %a = max(abs(u_m*(n_vec(1) + n_vec(2))), abs(u_p*(n_vec(1) + n_vec(2))));
    a=1;
    %a = abs(n_vec(1) + n_vec(2));
    %a = 0.1*abs(n_vec(1) + n_vec(2));
    F =    NaN(2,1);
    f_m = f(u_m); f_p = f(u_p);
    F(1) = (f_m+f_p);
    g_m = g(u_m); g_p = g(u_p);
    F(2) = (g_m+g_p);
    fhat = (1/2)*(dot(F,n_vec)) - a/2*(u_p-u_m);
end   

%%%%%%%%%%%%Start of LS functions%%%%%%%%%%%%%%%%%

function cell_avg = build_cell_avg_vec(inner_v_vec,f, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area)
    cell_avg = NaN(1,n2d);
    for ii = 1:length(inner_v_vec)
        i = inner_v_vec(ii);

        cell_int = int_over_cell(f,i,nod_in_e2d,tri,xcoord,ycoord);
        cell_avg(i) = cell_int/area(i);
    end
end

function  cell_int = int_over_cell(f,i,nod_in_e2d,tri,xcoord,ycoord)
    x_i = xcoord(i); y_i = ycoord(i);
    cell_int = 0;
    for j = 1:6
    big_t = nod_in_e2d(j,i);
    big_t_vert = tri(:,big_t);
    mid_points = NaN(2,2);
    count = 1;
    for v = 1:3
        if big_t_vert(v) ~= i
            mid_points(1,count) = (xcoord(big_t_vert(v)) + xcoord(i))/2;
            mid_points(2,count) = (ycoord(big_t_vert(v)) + ycoord(i))/2;
            count = count + 1;
        end
    end
    center = [sum(xcoord(big_t_vert))/3, sum(ycoord(big_t_vert))/3];
    
    lil_t_1 = [x_i, y_i; mid_points(:,1)'; center];
    lil_t_2 = [x_i, y_i; mid_points(:,2)'; center];
    
    cell_int = cell_int + tri_int(lil_t_1, f) +  + tri_int(lil_t_2, f);
    end
end

function result = tri_int(tri_coord, f)

%Updated to be symmetric, should work polynomials degree 2


    % tri_coord : 3 x 2 containing triangle coordinates to integrate over
    % f : function to integrate over

%Quad points on reference triangle

%xi_arr = [0, 1/2, 1/2];
xi_arr = [1/6,2/3,1/6];

%eta_arr = [1/2, 1/2, 0];
eta_arr = [1/6,1/6,2/3];

x_c = tri_coord(:,1); y_c = tri_coord(:,2);

%Compute area of triangle
A = 0.5*abs((x_c(1)*(y_c(2) - y_c(3)) + x_c(2)*(y_c(3) - y_c(1)) + x_c(3)*(y_c(1) - y_c(2))));

J = 2*A;

result = 0;

for i = 1:3
    xi = xi_arr(i); eta = eta_arr(i);

    x_p = (x_c(2) - x_c(1))*xi + (x_c(3) - x_c(1))*eta + x_c(1);

    y_p = (y_c(2) - y_c(1))*xi + (y_c(3) - y_c(1))*eta + y_c(1);

    result = result + f(x_p, y_p);
end

%Weight = 1/6
result = result *J/6;

end

function W = hexagon_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area)

    
        f1 = @(x,y) 1;    f2 = @(x,y) x;       f3 = @(x,y) y; 
        f4 = @(x,y) x.*y; f5 = @(x,y) x.^2; f6 = @(x,y) y.^2;
         %Build matrix for least squares
        A = zeros(m,k);
        for jj = 1:m %j is cell, jj is index
            j = nghbr_n2d(jj,i);
            area_j = area(j);
            A(jj,1) = 1/area_j*int_over_cell(f1,  j,nod_in_e2d,tri,xcoord,ycoord);
            A(jj,2) = 1/area_j*int_over_cell(f2,  j,nod_in_e2d,tri,xcoord,ycoord);
            A(jj,3) = 1/area_j*int_over_cell(f3,  j,nod_in_e2d,tri,xcoord,ycoord);
            A(jj,4) = 1/area_j*int_over_cell(f4,  j,nod_in_e2d,tri,xcoord,ycoord);
            A(jj,5) = 1/area_j*int_over_cell(f5,  j,nod_in_e2d,tri,xcoord,ycoord);
            A(jj,6) = 1/area_j*int_over_cell(f6  ,j,nod_in_e2d,tri,xcoord,ycoord);
        end
        
        %Build L matrix
        L = zeros(1,k);
        area_i = area(i);
        L(1) = 1/area_i*int_over_cell(f1,  i,nod_in_e2d,tri,xcoord,ycoord);
        L(2) = 1/area_i*int_over_cell(f2,  i,nod_in_e2d,tri,xcoord,ycoord);
        L(3) = 1/area_i*int_over_cell(f3,  i,nod_in_e2d,tri,xcoord,ycoord);
        L(4) = 1/area_i*int_over_cell(f4,  i,nod_in_e2d,tri,xcoord,ycoord);
        L(5) = 1/area_i*int_over_cell(f5,  i,nod_in_e2d,tri,xcoord,ycoord);
        L(6) = 1/area_i*int_over_cell(f6,  i,nod_in_e2d,tri,xcoord,ycoord);
        
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



 function [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d, edge_tri,edge_cross_dxdy] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d)
    % Function from code by Sergey Donilov
 n_num=zeros([1,n2d]);
 for n=1:e2d
    enodes=tri(:,n);
    n_num(enodes')=n_num(enodes')+1;
 end


 nod_in_e2d=zeros([max(n_num),n2d]);
 n_num=zeros([1,n2d]);
 for n=1:e2d   
    enodes=tri(:,n);
    n_num(enodes')=n_num(enodes')+1;
    for j=1:3
       nod_in_e2d(n_num(enodes(j)),enodes(j))=n;
    end
 end


  ng_num=zeros([1 n2d]); %Tracks num of neighbors in stencil
  check=zeros([1,n2d]);
  aux=zeros([1 10]);
  for j=1:n2d
     cc=0;
     for m=1:n_num(j)
         el=nod_in_e2d(m,j);
         for k=1:3
           a=tri(k,el);       
           if check(a)==0, 
              check(a)=1; 
              cc=cc+1;         
              aux(cc)=a;
           end
         end
	 
     end
     
     ng_num(j)=cc;
     check(aux(1:cc))=0;
  end


    % Repeat and fill   %%% TO DO %%%
  for j=1:n2d,
     cc=0;
     for m=1:n_num(j),
         el=nod_in_e2d(m,j);
         for k=1:3,
           a=tri(k,el);       
           if check(a)==0,  
              check(a)=1; 
              cc=cc+1;         
              aux(cc)=a;
           end;
         end;
     end;
     nghbr_n2d(1:cc,j)=aux(1:cc)';
     check(aux(1:cc))=0;      
  end;

  %Not sure if need all of this stuff yet... maybe for FV
[edges,edge_num, elem_edges,edge_tri,ed2d,ed2d_in,edge_dxdy, edge_cross_dxdy] = get_fv_arrays(tri, xcoord, ycoord ,n_num,ng_num, nghbr_n2d, n2d,nod_in_e2d,e2d);

   % rearrange neighbourhood so that n follows first %%% TO DO %%%
  for n=1:n2d,
     ng_num(n)=1;           % ng_num(n) is the number of neighbor vertices, including vertex n
     nghbr_n2d(1,n)=n;      % nghbr_n2d(1:ng_num(n),n) are their numbers.
  end;  

  for n=1:ed2d,
     n1=edges(1,n);
     n2=edges(2,n);
     nghbr_n2d(ng_num(n1)+1,n1)=n2;
     ng_num(n1)=ng_num(n1)+1;
     nghbr_n2d(ng_num(n2)+1,n2)=n1;
     ng_num(n2)=ng_num(n2)+1;
  end; 

% (a) Neighbourhood information 
 
 n_num=zeros([1,n2d]);
 for n=1:e2d
    enodes=tri(:,n);
    n_num(enodes')=n_num(enodes')+1;
 end;
 nod_in_e2d=zeros([max(n_num),n2d]);
 n_num=zeros([1,n2d]);
 for n=1:e2d   
    enodes=tri(:,n);
    n_num(enodes')=n_num(enodes')+1;
    for j=1:3,
       nod_in_e2d(n_num(enodes(j)),enodes(j))=n;
    end;
 end;

% (b) Neighbouring nodes:

  ng_num=zeros([1 n2d]);
  check=zeros([1,n2d]);
  aux=zeros([1 10]);
  for j=1:n2d,
     cc=0;
     for m=1:n_num(j),
         el=nod_in_e2d(m,j);
         for k=1:3,
           a=tri(k,el);       
           if check(a)==0,  
              check(a)=1; 
              cc=cc+1;         
              aux(cc)=a;
           end;
         end;
	 
     end;
     
     ng_num(j)=cc;
     check(aux(1:cc))=0;
  end;
  
  nghbr_n2d=zeros([max(ng_num),n2d]);   
  
  % Repeat and fill
  for j=1:n2d,
     cc=0;
     for m=1:n_num(j),
         el=nod_in_e2d(m,j);
         for k=1:3,
           a=tri(k,el);       
           if check(a)==0,  
              check(a)=1; 
              cc=cc+1;         
              aux(cc)=a;
           end;
         end;
     end;
     nghbr_n2d(1:cc,j)=aux(1:cc)';
     check(aux(1:cc))=0;      
  end;


  % ============
  % The areas of triangles:
  % ============
  elem_area=zeros([1,e2d]); 
  a=zeros([1,2]);
  b=zeros([1,2]);
  for n=1: e2d,
   elnodes=tri(:,n);
   a(1)=xcoord(elnodes(2))-xcoord(elnodes(1));
   a(2)=ycoord(elnodes(2))-ycoord(elnodes(1));
   b(1)=xcoord(elnodes(3))-xcoord(elnodes(1));
   b(2)=ycoord(elnodes(3))-ycoord(elnodes(1));
   elem_area(n)=0.5*abs(a(1)*b(2)-b(1)*a(2));
  end 
 
  % =============
  % Scalar cell (cluster) area in FV, 
  % coincides with lumped mass matriix in FE  
  % =============
 
  area=zeros([1,n2d]);
  for n=1:n2d,
    area(n)=sum(elem_area(nod_in_e2d(1:n_num(n),n)))/3.0;
  end 

 end



function order_plot(L2Error,Method,p,dx,err_name, k)
    figure(6)
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