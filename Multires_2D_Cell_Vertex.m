%Multi-resolution WENO implementation with cell-vertex mesh (constrained) 
    %%% Run pde = linear for order or pde = burgers for behavior around shock

clear all; close all;

pth = '/Users/chloegriffin/Brown Dropbox/Chloe Griffin/Research 2024/AWI Codes/Results';

    %Stored Info Guide: 
    % n2d: Number of nodes/control volumes/ triangle vertices (nx x ny)
    % e2d: Number of triangles

    % xcoord: (1xn2d) x-coordinates of control volume center/triangle vertices
    % ycoord  (1xn2d) y-coordinates of control volume center/triangle vertices
    
    % tri: (3xe2d) Tracks three vertices associated to each triangle (used
            % for plotting mesh)
    % nod_in_e2d: (6xn2d) Tracks triangles for each node (unordered)
    % nghbr_n2d:  (7xn2d) Tracks neighboring nodes for each node (unordered
            % for now) (could order)
    % ng_num: (1xn2d) Tracks num of neighbors for each node 
    % n_num:  (1xn2d) Tracks num of triagles for each node
    
    % edge_num/ed2d: Number of edges (not sure why both)
    % edges: (2xedge_num) Tracks vertices for each edge
    % edge_tri : (2xedge_num) Tracks triangle/triangles associated to each edge
    % elem_edges: (3xe2d) Tracks edges associated to each triangle (numbered as in
    %      edges)

    % edge_dxdy : (2xedge_num) Gives dx and dy for each edge
    % edge_cross_dxdy: (4xedge_num) Gives cross-distances for the edge
    % (used in gradient reconstruction)

    %enodes : tri(:,n) Nodes of triangle for a particular control volume

tic

%Uncomment the desired time
t_final = 2/2; %After shock
%t_final = 0.3/2; 

%Choose eps to test
eps = 1e-5;
%t_final = 1e-2;
avoid_shock = true;

p = 3;


approach = "standard";
%approach = "beta_1_0";
%approach = "single_beta_ls"; 
%approach = "dir_deriv_sten_1";
%approach = "dir_deriv_sten_2";

pde = "burgers";
%pde = "linear";
alpha = 1/3; beta = 2/3; 
u_initial = @(x,y) alpha + beta*sin((x+y)); 

c=0.5;

if pde == "burgers"
f = @(u) u.^2/2;
g = @(u) u.^2/2;
elseif pde == "linear"
f = @(u) c*u;
g = @(u) c*u; 
end

  n_array = [20,40,80];

   Method = "Multiresolution WENO"; p = 3;
        

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
          n_x = n_set_array(kk); 

          h=(xR-xL)/n_x;
          x = (xL:h:xR);
          y = (yL:h:yR); 
        
          dx  = x(2:end) - x(1:end-1);
          dy  = y(2:end) - y(1:end-1);
        
          %Will update later for nonuniform case
          
          x = [x(1) - dx(end) - dx(end - 1) - h,x(1) - dx(end) - dx(end - 1) , x(1) - dx(end), x, x(end) + h, x(end) + dx(1) + dx(2), x(end) + dx(1) + dx(2)+h];
          y = [y(1) - dy(end) - dy(end - 1) - h,y(1) - dy(end) - dy(end - 1),  y(1) - dx(end), y, y(end) + h, y(end) + dy(1) + dy(2), y(end) + dy(1) + dy(2)+h];

          N_x = length(x);
          N_y = length(y);

          % n_x = N_x - 7;
          % n_y = N_y - 7;

         DX(kk)=(xR-xL)/(n_x);
       
         [tri, nodnum, xcoord, ycoord, e2d, n2d] = gen_reg_tri_mesh(x,y);

         [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d, edge_tri,elem_edges,edge_cross_dxdy] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d);

          % area_in = reshape(area, N_y,N_x);
          % area_in = area_in(4:end - 3, 4:end - 3);
          % area_in = reshape(area_in, 1,(N_y-6)*(N_x-6));


         [inner_v_vec_2,inner_v_vec,ind_vec,n_steps,inner_steps] = get_indexing_arrays(N_x,N_y);

         nghbr_n2d = order_neighbor(nghbr_n2d,nod_in_e2d, edges,edge_tri,tri,elem_edges,inner_v_vec_2, inner_steps);

         cell_avg = build_cell_avg_vec(inner_v_vec, u_initial, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);

         cell_avg_array = reshape(cell_avg, N_y, N_x);

        %Third order approximation
         k = 6; m = 7;
        

          W = NaN(k, m+1, N_x*N_y);
          for ii = 1:inner_steps
            i = inner_v_vec_2(ii);
            W(:,:,i) = cell_vertex_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area);
          end

          B = NaN(m+1, m+1, N_x*N_y);
          for ii = 1:inner_steps
            i = inner_v_vec_2(ii);
             B(:,:,i) = cell_vertex_smoothness_matrix(squeeze(W(:,:,i)),i,nod_in_e2d,tri,xcoord,ycoord,area);
          end

          L1 = NaN(m+1, m+1, N_x*N_y);
          for ii = 1:inner_steps
            i = inner_v_vec_2(ii);
             L1(:,:,i) = lin_cell_vertex_smoothness_matrix(squeeze(W(:,:,i)),i,nod_in_e2d,tri,xcoord,ycoord,area);
          end

          %Precompute quad points

          edge_quad = NaN(4,2,ed2d);
          for e = 1:ed2d
            v1 = edges(1,e); v2 = edges(2,e);

            %Compute results for interior points (where two triangles adjacent) 
            %Alternatively can run until edge_in

            if (edge_tri(1,e) && edge_tri(2,e) ~= 0)
                cross_info = edge_cross_dxdy(:,e);
     
                midpoint = [(xcoord(v1) + xcoord(v2))/2, (ycoord(v1) + ycoord(v2))/2]; %Compute midpoint
        
                for side = 1:2
                    cross_vec = cross_info((2*side - 1): 2*side)';                
                    edge_quad((2*side - 1): 2*side,1,e) = midpoint + cross_vec*(3 - sqrt(3))/6; edge_quad((2*side - 1): 2*side,2,e) = midpoint + cross_vec*(3 + sqrt(3))/6;            
                end
            end
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
            u_1 = u + dt*L(u,n_steps, N_x, N_y, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,B,a,edge_cross_dxdy,eps,approach,edge_quad,L1,pde);
            u_1 = apply_bc(u_1,N_x,N_y);

            u_2 = 3/4*u +1/4*u_1 + 1/4*dt*L(u_1,n_steps, N_x,  N_y, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,B,a,edge_cross_dxdy,eps,approach,edge_quad,L1,pde);
            u_2 = apply_bc(u_2,N_x,N_y);

            u_next = 1/3*u + 2/3*u_2 + 2/3*dt*L(u_2,n_steps, N_x,  N_y, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,B,a,edge_cross_dxdy,eps,approach,edge_quad,L1,pde);
            u_next = apply_bc(u_next,N_x,N_y);
    
            u = u_next;
            t = t + dt;
        end
        toc
        %Prep for visualization

        if pde == "burgers"
            u_ref= @(x,y) solveBurgers2(alpha, beta, 2*t_final,(x+y));

         elseif pde == "linear"
            u_ref = @(x,y) u_initial(x - c*t_final, y - c*t_final);
        end

        x_dense = (xL:h/10:xR)';
        y_dense = (yL:h/10:yR);
        u_ref_x = u_ref(x_dense, 0.001);
        u_ref_y = u_ref(0.001, y_dense);
        fig_num = kk;
        fig_num = fig_num+1; variable_name = "sin(x+y)"; 
        Method = "Least Squares Poly Approx"; 

        if n_x == n_set_array(end)
             %Plot exact cell averages
             exact_avg_vec = build_cell_avg_vec(inner_v_vec, u_ref, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
             fig_num = 3;
             exact_avg_array = reshape(exact_avg_vec,N_y,N_x);
             figure(fig_num)
             trimesh(tri', xcoord, ycoord, exact_avg_array)
             title("Plot of Exact Cell Avg at t=" + t_final + " with nX = "+ n_x + " and eps = " + eps)
             saveas(fig_num, sprintf('%s/2D/hexagon_multres_tFinal_%d_eps_%d_plot_%s_exact.jpg',pth,t_final*2,eps,approach),'jpeg')

             view([-1, 1, 0]);
             saveas(fig_num, sprintf('%s/2D/hexagon_multres_tFinal_%d_eps_%d_plot_exact_%s_side.jpg',pth,t_final*2,eps,approach),'jpeg')

             %Plot approx cell averages
             fig_num = 5;
             approx_avg_array = reshape(u,N_y,N_x);
             figure(fig_num)
             trimesh(tri', xcoord, ycoord, approx_avg_array)
             title("Plot of Approx Cell Avg at t=" + t_final + " with nX = " + n_x + " and eps = " + eps)
             saveas(fig_num, sprintf('%s/2D/hexagon_multres_tFinal_%d_eps_%d_plot_approx.jpg',pth,t_final*2,eps),'jpeg')
             saveas(fig_num, sprintf('%s/2D/hexagon_multres_tFinal_%d_eps_%d_plot_approx_%s.jpg',pth,t_final*2,eps,approach),'jpeg')

             view([-1, 1, 0]);
             saveas(fig_num, sprintf('%s/2D/hexagon_multres_tFinal_%d_eps_%d_plot_approx_%s_side.jpg',pth,t_final*2,eps,approach),'jpeg')

             fig_num = fig_num + 3;
        end
    
        exact_avg_vec = build_cell_avg_vec(inner_v_vec, u_ref, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
       
        exact_avg_array = reshape(exact_avg_vec,N_y,N_x);
        approx_avg_array = reshape(u,N_y,N_x);

        if avoid_shock == false
        uDif = approx_avg_array(4:end - 3, 4:end -3) - exact_avg_array(4:end - 3, 4:end - 3);

        %L2_Error_quadratic(kk) = ((1/(n_steps))*sum(sum(abs(uDif).^2)))^(1/2);
        
        elem_area = reshape(area,N_y,N_x);
        elem_area = elem_area(4:end-3,4:end-3);
        L2_Error_quadratic(kk) = (sum(sum(elem_area.*abs(uDif).^2)))^(1/2);

        end

        if avoid_shock == true

        %Find shock idx (which cells are outside of shock)

            global_L2_Error = 0;
  
            xShock = 2*pi + alpha * t_final;
            diag_shock1 = xShock - pi;
            diag_shock2 = xShock + pi;    
            
            dist_shock1 = abs(xcoord + ycoord - diag_shock1);
            dist_shock2 = abs(xcoord + ycoord - diag_shock2);

            padding = 1;

            in_shock1 = dist_shock1 <= padding;
            in_shock2 = dist_shock2 <= padding;

            away_from_shock = ~(in_shock1 |  in_shock2);
            cell_count = sum(away_from_shock);
           
            away_from_shock = reshape(away_from_shock, N_y,N_x);

            exact_avg_array = exact_avg_array(4:end-3, 4:end-3);
            approx_avg_array = approx_avg_array(4:end-3, 4:end-3);
            away_from_shock = away_from_shock(4:end-3, 4:end-3);

            uDif = approx_avg_array(away_from_shock) - exact_avg_array(away_from_shock);

            L2_Error_quadratic(kk) = ((1/(cell_count))*sum(sum(abs(uDif).^2)))^(1/2);

        end

   
   end
   
    format shorte
    p = 3;
    err_name = "L2 Error";
    T1 = order_table(DX,n_array, L2_Error_quadratic, err_name);
    order_plot(L2_Error_quadratic, Method,p,DX,err_name, k);
    saveas(10, sprintf('%s/2D/hexagon_multres_tFinal_%d__eps_%d_convergence_plot_%s.jpg',pth,t_final*2,eps,approach),'jpeg')


%%%%%%%%%%%Start of new functions%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inner_v_vec_2,inner_v_vec,ind_vec,n_steps,inner_steps] = get_indexing_arrays(N_x,N_y)

 %Useful in implementing periodic bounday condition or using ghost points
 % N_x = length(x);
 % N_y = length(y);

 %Provides 3 boundary layers: 
 %  inner_v_array: ensures control volumes have full mass
 %  inner_v_array_2: ensures that we have reconstruction with full stencil
 %  ind_array: Has reconstructed neighbors for flux computation


 bound_tri_vec = 1:N_x*N_y;
 bound_tri_array = reshape(bound_tri_vec, N_y,N_x);

 inner_v_array = bound_tri_array(2:end - 1, 2:end - 1); %Array for all cell averages
 inner_v_array_2 = inner_v_array(2:end-1, 2:end-1);     %Array for all reconstruction
 ind_array =       inner_v_array(3:end - 2, 3:end - 2); %Array for stored computation

 inner_v_vec = reshape(inner_v_array, (N_x-2)*(N_y-2),1);    %Vector for all cell averages
 inner_v_vec_2 = reshape(inner_v_array_2, (N_x-4)*(N_y-4),1);%Vector for all reconstruction
 ind_vec = reshape(ind_array, (N_x - 6)*(N_y - 6),1);        %Vector for stored computation

  n_steps = (N_x - 6)*(N_y - 6); %Size of ind_vec
  inner_steps = (N_x-4)*(N_y-4); %Size of inner_v_vec_2

end


function ordered_nghbr_n2d  = order_neighbor(nghbr_n2d,nod_in_e2d, edges,edge_tri,tri,elem_edges,inner_v_vec_2, inner_steps)
    ordered_nghbr_n2d = nghbr_n2d;

    for jj = 1:inner_steps
        v_j = inner_v_vec_2(jj);

        %Set first entry as v_j
        ordered_nghbr_n2d(1,v_j) =  v_j; 

        tri_1 = nod_in_e2d(1,v_j);
        %Get edges attached to vertex

        e_tri_1 = elem_edges(:,tri_1);
        v_tri_1=  edges(:,e_tri_1);

        % Check presence of v_j in each edge
        contains_x = [any(v_tri_1(:,1) == v_j), any(v_tri_1(:,2)  == v_j), any(v_tri_1(:,3)  == v_j)];
       
        % Find vertices of edges that contain v_j
        indices =  find(contains_x);
        v_edge_1 = edges(:,e_tri_1(indices(1)));
        v_edge_2 = edges(:,e_tri_1(indices(2)));

        e_2 = e_tri_1(indices(2));

        %Assign correct vertices in order 
        if v_edge_1(1) == v_j
            ordered_nghbr_n2d(2,v_j) =  v_edge_1(2);
        else
            ordered_nghbr_n2d(2,v_j) =  v_edge_1(1);
        end


        %Assign correct vertices in order 
        if v_edge_2(1) == v_j
            ordered_nghbr_n2d(3,v_j) =  v_edge_2(2);
        else
            ordered_nghbr_n2d(3,v_j) =  v_edge_2(1);
        end

        e_m = e_2;
        tri_mm1 = tri_1;

        for m = 2:5
            check_tri_m = edge_tri(:,e_m);
            if check_tri_m(1) == tri_mm1
                tri_m = check_tri_m(2);
            else
                tri_m = check_tri_m(1);
            end
        
            e_tri_m = elem_edges(:,tri_m);
            v_tri_m=  edges(:,e_tri_m);
    
            % Check presence of v_j in each edge
            contains_x = [any(v_tri_m(:,1) == v_j), any(v_tri_m(:,2)  == v_j), any(v_tri_m(:,3)  == v_j)];
           
            % Find vertices of edges that contain v_j
            indices =  find(contains_x);
    
            check_e_mp1 = e_tri_m(indices);
    
            if check_e_mp1(1) == e_m
                    e_mp1= check_e_mp1(2);
            else
                    e_mp1= check_e_mp1(1);
            end

            v_edge_mp1 = edges(:,e_mp1);

            if v_edge_mp1(1) == v_j
                ordered_nghbr_n2d(m+2,v_j) =  v_edge_mp1(2);
            else
                ordered_nghbr_n2d(m+2,v_j) =  v_edge_mp1(1);
            end


            e_m = e_mp1;
            tri_mm1 = tri_m;
        end


    end



end



function u_1 = apply_bc(u_1,N_x,N_y)
u_1 = reshape(u_1,N_y,N_x);
%Pad twice periodically

u_1(3,:) = u_1(end - 4,:);
u_1(end-2,:) = u_1(5,:);
u_1(:,3) = u_1(:,end - 4);
u_1(:,end-2) = u_1(:,5);

u_1(2,:) = u_1(end - 5,:);
u_1(end-1,:) = u_1(6,:);
u_1(:,2) = u_1(:,end - 5);
u_1(:,end-1) = u_1(:,6);
u_1 = reshape(u_1,1,N_x*N_y);
end


function B = cell_vertex_smoothness_matrix(W,i, nod_in_e2d,tri,xcoord,ycoord,area)
   
    f1 = @(x,y) 1;    f2 = @(x,y) x;       f3 = @(x,y) y; 
    f4 = @(x,y) x.*y; f5 = @(x,y) x.^2; f6 = @(x,y) y.^2;

    L1 = int_over_cell(f1,i,nod_in_e2d,tri,xcoord,ycoord);
    L2 = int_over_cell(f2,i,nod_in_e2d,tri,xcoord,ycoord);
    L3 = int_over_cell(f3,i,nod_in_e2d,tri,xcoord,ycoord);
    L4 = int_over_cell(f4,i,nod_in_e2d,tri,xcoord,ycoord);
    L5 = int_over_cell(f5,i,nod_in_e2d,tri,xcoord,ycoord);
    L6 = int_over_cell(f6,i,nod_in_e2d,tri,xcoord,ycoord);

    W2 = W(2,:); W3 = W(3,:); 
    W4 = W(4,:); W5 = W(5,:); W6 = W(6,:);

    area_i = area(i);

    B1 = ((W2')*W2 + (W3')*W3 + area_i*(W4'*W4 + 4*(W5')*W5 + 4*(W6')*W6))*L1;
    B2 = (4*(W2')*W5 + 2*(W3')*W4)*L2;
    B3 = (2*(W2')*W4 + 4*(W3')*W6)*L3;
    B4 = (4*(W4')*W5 + 4*(W4')*W6)*L4;
    B5 = (4*(W5')*W5 + (W4'*W4))*L5;
    B6 = ((W4')*W4 + 4*(W6')*W6)*L6;

    %B = 121/100*(B1 + B2 + B3 + B4 + B5 + B6) ;
    B = (B1 + B2 + B3 + B4 + B5 + B6);

end

function L1 = lin_cell_vertex_smoothness_matrix(W,i, nod_in_e2d,tri,xcoord,ycoord,area)
   
    W2 = W(2,:); W3 = W(3,:); 

    area_i = area(i);

    L1 = ((W2')*W2 + (W3')*W3)*area_i;
end


function u_hold = L(u,n_steps, N_x,N_y, ed2d,xcoord,ycoord,edges,edge_tri,tri,W,f,g,nghbr_n2d,ind_vec,area,inner_v_vec_2,inner_steps,B,a,edge_cross_dxdy,eps,approach,edge_quad,L1,pde)
    %quad_res = weno_reconstruction(W,sten_array,coord_array,quad_vals,v_vec,N); 
    %flux_array = find_fluxes(quad_res,sten_array,N,f,g,a);

    inner_sum = zeros(1,N_x*N_y);

    %For LS reconstruction, store coefficients and weights for all i 
    P = weno_reconstruction(u,W,nghbr_n2d,N_x,N_y,inner_v_vec_2,inner_steps,B,xcoord,ycoord,area,eps,approach,L1);
    
    u_cell = @(x,y,i) P(1,i)*u(i)+ P(2,i)*(11/10*(P(3,i) + P(4,i)*x + P(5,i)*y+ P(6,i)*x.*y + P(7,i)*x.^2 + P(8,i)*y.^2) - 1/10*u(i));
    
    %Uncomment to test fully quadratic or linear
        %u_cell = @(x,y,i) P(3,i) + P(4,i)*x + P(5,i)*y+ P(6,i)*x.*y + P(7,i)*x.^2 + P(8,i)*y.^2;
        %u_cell = @(x,y,i) u(i);
    
    for e = 1:ed2d
        v1 = edges(1,e); v2 = edges(2,e);
    
        %Compute results for interior points (where two triangles adjacent) 
            %Alternatively can run until edge_in
    
        if (edge_tri(1,e) && edge_tri(2,e) ~= 0)
            cross_info = edge_cross_dxdy(:,e);
    
            for side = 1:2
                cross_vec = cross_info((2*side - 1): 2*side)';
                quad_points = edge_quad((2*side - 1):2*side,:,e);
    
                %Get outward facing normal vector to side edge
                %   Rotate 90 degrees/-90 degrees
                if side == 1
                     n_vec = [cross_vec(2); -cross_vec(1)];
    
                else 
                     n_vec = [-cross_vec(2); cross_vec(1)];
    
                end
    
                %Vectorized computation with quad points
                 
                %Use LS/multi-res for scalar reconstruction from either cell
                uq_m = u_cell(quad_points(1,:),quad_points(2,:),v1);
                uq_p = u_cell(quad_points(1,:),quad_points(2,:),v2);
            
                %fv_1(q) = lf(uq_1,uq_2,n_vec,f,g); %u_m, u_p
    
                %LF
                    %burgers:
                    if pde == "burgers"
                    a = max(abs(uq_m.*(n_vec(1) + n_vec(2))), abs(uq_p.*(n_vec(1) + n_vec(2))));
                    elseif pde == "linear"
                    a = ones(1,2)*abs(n_vec(1) + n_vec(2));
                    end

                    F_m = [f(uq_m); g(uq_m)]; %2x2
                    F_p = [f(uq_p); g(uq_p)]; %2x2
                    fv_1 = n_vec'*(F_m + F_p) - a.*(uq_p-uq_m); %1/2 at end
                    fv_2 = -fv_1;
          
    
                inner_sum(v1) = inner_sum(v1) + sum(fv_1); %Summing over quad points for edge
                inner_sum(v2) = inner_sum(v2) + sum(fv_2); %weights = 1/2
             end
         end
       
    end

    % 
    % u_hold = NaN(1,N_x*N_y);
    % for ii = 1:n_steps
    %     i = ind_vec(ii);
    %     area_i = area(i);
    %     u_hold(i) = -(1/area_i)*inner_sum(i);
    % end

    u_hold = -inner_sum./(4*area);
end


function P  = weno_reconstruction(u,W,nghbr_n2d,N_x,N_y,inner_v_vec_2,inner_steps,B,xcoord,ycoord,area,eps,approach,L1)
      
      %Computing polynomial weights and coefficients for reconstuction in
      %all cells as indexed by inner_v_vec_2

      k = 6; m = 7;
      P =  NaN(k+2,N_x*N_y);


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
    
        C = W(:,:,i)*d;

        %q1 = @(x,y) v_vec(i);
        %q2 = @(x,y) C(1) + C(2)*x + C(3)*y+ C(4)*x.*y + C(5)*x.^2 + C(6)*y.^2;
    
        %p1 = @(x,y) q1(x,y);
        %p2 = @(x,y) (11/10)*q2(x,y) - (1/10)*p1(x,y);
    
        %Build B_1 (Linear smoothness from three different stencils)
        area_i = area(i);


        stencil_1 = [1,3,6];
        stencil_2 = [1,7,4];
        stencil_3 = [1,2,5];


        if approach == "standard" || approach == "dir_deriv_sten_2"
        stencil_1 = [1,3,4];
        stencil_2 = [1,5,6];
        stencil_3 = [1,7,2];
        end

        if approach ~= "beta_1_0" && approach ~= "single_beta_ls"

        [C2_1, C3_1, b1_1] = compute_lin_smooth(i,stencil_1,u,nghbr_n2d,xcoord,ycoord,area_i,approach);
        [C2_2, C3_2, b1_2] = compute_lin_smooth(i,stencil_2,u,nghbr_n2d,xcoord,ycoord,area_i,approach);
        [C2_3, C3_3, b1_3] = compute_lin_smooth(i,stencil_3,u,nghbr_n2d,xcoord,ycoord,area_i,approach);


        beta_dif_avg =(abs(b1_1 - b1_2)+abs(b1_1 - b1_3) + abs(b1_2 - b1_3))/3;

        sigma_1 = 1/3*(1 + (beta_dif_avg)^2/(b1_1 + eps));
        sigma_2 = 1/3*(1 + (beta_dif_avg)^2/(b1_2 + eps));
        sigma_3 = 1/3*(1 + (beta_dif_avg)^2/(b1_3 + eps));

        sigma = sigma_1 + sigma_2  + sigma_3;


        beta_1 = (area_i/sigma^2)*(sigma_1*(C2_1) + sigma_2*(C2_2)+...
                               + sigma_3*(C2_3))^2+...
                +(area_i/sigma^2)*(sigma_1*(C3_1) + sigma_2*(C3_2)+...
                               + sigma_3*(C3_3))^2;

         %Alternative to linear combo:
             %beta_1 = min([b1_1,b1_2,b1_3]); 


        elseif approach == "beta_1_0"
            beta_1 = 0;
        elseif approach == "single_beta_ls"
            beta_1 = d'*L1(:,:,i)*d;
        end

      

        gam_12 = 1/11; gam_22 = 10/11; 

        %Build B_2 using precomputed matrix
            %Note: Using smoothness indicator from q only
        beta_2 = d'*B(:,:,i)*d;
   

        %Compute non-linear weights (WENO-Z)
        tau_2 = abs(beta_2 - beta_1);
   
        omega_bar_12 = gam_12*(1 + tau_2/(eps + beta_1));
        omega_bar_22 = gam_22*(1 + tau_2/(eps + beta_2));
    
        omega_12 = omega_bar_12/(omega_bar_12 + omega_bar_22);
        omega_22 = omega_bar_22/(omega_bar_12 + omega_bar_22);

        P(1,i) = omega_12; P(2,i) = omega_22; P(3:8,i) = C(:);
      end 
end


function [C2, C3, b1] = compute_lin_smooth(i,stencil,u,nghbr_n2d,xcoord,ycoord,area_i,approach)
        eps = 1e-12;

        cell_1 = i; cell_2 = nghbr_n2d(stencil(2),i); cell_3 = nghbr_n2d(stencil(3),i);

        v1 = u(cell_1); v2 = u(cell_2); v3 = u(cell_3);

        %Note: for non-uniform, center vertex != centroids; approach =
        %standard could break down

        xc_1 = xcoord(cell_1); xc_2 = xcoord(cell_2); xc_3 = xcoord(cell_3);
        yc_1 = ycoord(cell_1); yc_2 = ycoord(cell_2); yc_3 = ycoord(cell_3);

        dx2 = xc_2 - xc_1 ; dx3 = xc_3 - xc_1;
        dy2 = yc_2 - yc_1;  dy3 = yc_3 - yc_1;
        dv2 = v2 - v1;      dv3 =  v3 - v1;

        %Note: Could check if determinant is zero to see if colinear
        
            if approach == "dir_deriv_sten_1" ||  approach == "dir_deriv_sten_2"

            C2 = 0; C3 = 0;
            C_arr = [0,0];
            check = 0;
            if abs(dx2) > eps
                C_arr(1) = dv2/dx2;
                check = check + 1;
            elseif abs(dx3) > eps
                C_arr(2) = dv3/dx3;
                check = check + 1;
            end

            if check == 2
                C2 = min(C_arr);
                %C2 = mean(C_arr);
            elseif check == 1
                C2 = max(C_arr);
            end


            C_arr = [0,0];
            check = 0;
            if abs(dy2) > eps
                C_arr(1) = dv2/dy2;
                check = check + 1;
            elseif abs(dy3) > eps
                C_arr(2) = dv3/dy3;
                check = check + 1;
            end

            if check == 2
                C3 = min(C_arr);
                %C3 = mean(C_arr);
            elseif check == 1
                C3 = max(C_arr);
            end
            end


        %Standard Uniform Approach
        if approach == "standard"
        rhs = [dv2; dv3];

        A = [dx2, dy2; dx3, dy3];

        C = A\rhs;

        C2 = C(1); C3 = C(2);
        end

        %Add check 
        % P = @(x,y)C(1)*x+C(2)*y;
        % if (P(xc_2, yc_2) - P(xc_3, yc_3)) - (v2 - v3) > eps
        %     error("Linear poly doesn't fit")
        % end



        %Full solve (doesn't work if any colinear)
        % A = [1, xc_1, yc_1; 1, xc_2, yc_2; 1, xc_3, yc_3];
        % v_vec = [v1;v2;v3];
        % C = A\v_vec;
        % 
        % C2 = C(2); C3 = C(3);

        b1 = (C2^2+C3^2)*area_i;
end

function fhat = lf(u_m,u_p,n_vec,f,g,a)
    %burgers:
    a = max(abs(u_m*(n_vec(1) + n_vec(2))), abs(u_p*(n_vec(1) + n_vec(2))));
        %a = max(abs(u_m), abs(u_p));

    %linear:
        %a=1;

    % F =    NaN(2,1);
    % f_m = f(u_m); f_p = f(u_p);
    % F(1) = (f_m+f_p);
    % g_m = g(u_m); g_p = g(u_p);
    % F(2) = (g_m+g_p);
    % fhat = (1/2)*(dot(F,n_vec) - a*(u_p-u_m));


    F_m = [f(u_m); g(u_m)];
    F_p = [f(u_p); g(u_p)];
    fhat = (1/2)*(dot((F_m + F_p),n_vec) - a*(u_p-u_m));
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

%Approach from literature, only accurate for monomials, not polynomials,
%  Updated to be symmetric, should work polynomials degree 2

    %%%tri_coord : 3 x 2 containing triangle coordinates to integrate over


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
    
    result = result *J/6;

end

function W = cell_vertex_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area)

    
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



 function [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d, edge_tri,elem_edges,edge_cross_dxdy] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d)

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
    figure(10)
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