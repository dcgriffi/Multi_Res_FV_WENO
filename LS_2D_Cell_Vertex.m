%Convergence test for LS 3rd order polynomials with uniform cell-vertex mesh (scaled, constrained) 


clear all; close all;
format shorte;

    %Stored Info Guide: 
    % n2d: Number of nodes/control volumes/ triangle vertices (nx x ny)
    % e2d: Number of triangles

    % xcoord: (1xn2d) x-coordinates of control volume center/triangle vertices
    % ycoord  (1xn2d) y-coordinates of control volume center/triangle vertices
    
    % tri: (3xe2d) Tracks three vertices/control volumes associated to each triangle (used
            % for plotting mesh)
    % nod_in_e2d: (6xn2d) Tracks triangles for each node (unordered)
    % nghbr_n2d:  (7xn2d) Tracks neighboring nodes for each node (unordered
            % for now)
    % ng_num: (1xn2d) Tracks num of neighbors for each node 
    % n_num:  (1xn2d) Tracks num of triagles for each node
    
    % edge_num/ed2d: Number of edges (not sure why both) (Not sure if
             % triangle or control volume edges now)
    % edges: (2xedge_num) Tracks vertices for each edge
    % edge_tri : (2xedge_num) Tracks triangle/triangles associated to each edge
    % elem_edges: (3xe2d) Tracks edges associated to each triangle (numbered as in
    %      edges)

    % edge_dxdy : (2xedge_num) Gives dx and dy for each edge (from v1 to v2)
    % edge_cross_dxdy: (4xedge_num) Gives cross-distances for the edge
    % (from midpoint to centroid of triangle)

    %enodes : tri(:,n) Nodes of triangle for a particular control volume
 
   %n_array = [4,8,16,32,64,128];
   %n_array = [5,9,17,33,65,129,257];
   %n_array = [33,65,129,257];
   %n_array = [11,21,41,81,161];
   %n_array = [129,257];

   n_array = [20,40,80,160];

   Method = "Least Squares Poly Approx"; p = 3;
        

  %Preallocate arrays for memory
    DX = zeros(length(n_array),1); dt = zeros(length(n_array),1); 
    L2_Error_quadratic = zeros(length(n_array),1); L1Error = zeros(length(n_array),1); 
    L2_Error_linear = zeros(length(n_array),1);
    LinfError= zeros(length(n_array),1); err = zeros(length(n_array),1); 

  n_set_array = n_array;

  xL = 0; xR = 2*pi;
  yL = 0; yR = 2*pi;


  for kk = 1:length(n_array)
  n = n_set_array(kk); 
  h=(xR-xL)/n;
  x = (xL:h:xR);
  y = (yL:h:yR); 

  dx  = x(2:end) - x(1:end-1);
  dy  = y(2:end) - y(1:end-1);

  %Will update later for nonuniform case
  % x = [x(1) - 2*h, x(1) - h, x, x(end) + h, x(end) + 2*h];
  % y = [y(1) - 2*h, y(1) - h, y, y(end) + h, y(end) + 2*h ];

  x = [x(1) - dx(end) - dx(end - 1), x(1) - dx(end), x, x(end) + h, x(end) + dx(1) + dx(2)];
  y = [y(1) - dy(end) - dy(end - 1), y(1) - dx(end), y, y(end) + h, y(end) + dy(1) + dy(2)];
        
 N = n + 5;
 DX(kk)=(xR-xL)/N;
% 
% nx = length(x)-1;
% ny = length(y)-1;
% N = (nx)*(ny);

 alpha = 1/3; beta = 2/3;
 u_initial = @(x,y) alpha + beta*sin((x+y)/2); 


 % N = 8;
 % 
  %x = 1:1:8;
  %y = 1:1:8;

 [tri, nodnum, xcoord, ycoord, e2d, n2d] = gen_reg_tri_mesh(x,y);

 [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d);

 bound_tri_vec = 1:N*N;

 bound_tri_array = reshape(bound_tri_vec, N,N);
 inner_v_array = bound_tri_array(2:end - 1, 2:end - 1);
 ind_array = inner_v_array(2:end - 1, 2:end - 1);

 inner_v_vec = reshape(inner_v_array, (N-2)*(N-2),1);
 ind_vec = reshape(ind_array, (N-4)*(N-4),1);

 n_steps = (N-4)*(N-4);

 %Build inner ind_array (get reconstruction of these)

 cell_avg = build_cell_avg_vec(inner_v_vec, u_initial, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
 
 cell_avg_array = reshape(cell_avg, N, N);

 v_array = cell_avg_array(3:end - 2, 3:end - 2);
 v_vec = reshape(v_array, N-4,N-4); %This is what I would actually want to plot and update in time stepping

 if n+1== n_array(2)
 %Plot just mesh
 figure(1)
trimesh(tri', xcoord', ycoord', zeros(N*N,1),"EdgeColor","c"); 
view(2)
title("Plot of Mesh with nX = " + N)

 %Plot cell_avg_array
 figure(2)
 trimesh(tri', xcoord, ycoord, cell_avg_array)
 title("Plot of Cell Averages with nX = " + N)
 end

 k = 6; m = 7;

 W = zeros(k, m+1, n_steps);
 for ii = 1:n_steps
   i = ind_vec(ii);
   W(:,:,ii) = hexagon_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area);
 end

  C = zeros(n_steps,k);
  %Build coefficient matrix using weight matrix
  for ii = 1:n_steps
    i = ind_vec(ii);
    stencil = nghbr_n2d(:,i);
    d = zeros(m+1,1);
    for j = 1:m
        cell = stencil(j);
        d(j) = cell_avg(cell);
    end
    d(m+1) = -cell_avg(i);

    C(ii,:) = W(:,:,ii)*d;
  end


  ref_avg = build_ref_array(C,area,n_steps, ind_vec,nod_in_e2d,tri,xcoord,ycoord,n2d);
     
  ref_avg_array = reshape(ref_avg,N,N);



  global_L2_quadratic = global_L2_error(C,u_initial,n_steps,ind_vec,nod_in_e2d,tri,xcoord,ycoord);

  %uDif = cell_avg_array(4:end - 3, 4:end -3) - ref_avg_array(4:end - 3, 4:end - 3);

  % global_L2_quadratic = (1/(n_steps))*sum(sum(abs(uDif)).^2);
  % 
  % %global_L2_quadratic = 
   L2_Error_quadratic(kk)= (global_L2_quadratic)^(1/2);

  %Add visualization of least square result
   
     if n+1== n_array(end - 1)
     ref_avg = build_ref_array(C,area,n_steps, ind_vec,nod_in_e2d,tri,xcoord,ycoord,n2d);
     
     ref_avg_array = reshape(ref_avg,N,N);
     figure(3)
     trimesh(tri', xcoord, ycoord, ref_avg_array)
     title("Plot of Least Square Cell Averages with nX = " + N)

     check_array = abs(cell_avg_array - ref_avg_array);
     check_array  = check_array(3:end - 2, 3:end - 2);
     check = max(max(check_array))
     end





     %Add second visualization of least square result: Plot the error on
     %the grid, for a single C value

     if n+1== n_array(1)

     %Cell of interest = cell_i

    % cell_i = 28; %second from bottom and second from left %Qu:How to map from full array to ind_vec
    % cell_i = 29; %second from bottom and second from left %Qu:How to map from full array to ind_vec

      %ii = find(ind_vec==cell_i);

      ii = 6;

      f = @(x,y) C(ii,1) + C(ii ,2)*x + C(ii ,3)*y+ C(ii ,4)*x.*y + C(ii ,5)*x.^2 + C(ii ,6)*y.^2;
        
      fix_avg = build_cell_avg_vec(inner_v_vec, f, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
     
      fix_avg_array = reshape(fix_avg,N,N);
      figure(4)
      trimesh(tri', xcoord, ycoord, fix_avg_array)
      title("Plot of Fixed Poly at " + ii + " with nX = " + N)

      fix_avg_array = reshape(fix_avg,N,N);
      figure(5)
      trimesh(tri', xcoord, ycoord, cell_avg_array)
      title("Plot of Correct Cell Avg at " + ii + " with nX = " + N)

      check_array = abs(cell_avg_array - fix_avg_array);
      check_array  = check_array(3:end - 2, 3:end - 2)
      %check = max(max(check_array))
      end
        


  end

  
  err_name = "L2 Error";
  T1 = order_table(DX,n_array, L2_Error_quadratic, err_name);
  order_plot(L2_Error_quadratic, Method,p,DX,err_name, k);


function global_L2_Error = global_L2_error(C,exact_sol,n_steps,ind_vec,nod_in_e2d,tri,xcoord,ycoord)
    global_L2_Error = 0;
    for jj = 1:n_steps
        cell = ind_vec(jj);
        %Obtain approx function at each cell (change if increase k)
        
        u_cell = @(x,y) C(jj ,1) + C(jj ,2)*x + C(jj ,3)*y+ C(jj ,4)*x.*y + C(jj ,5)*x.^2 + C(jj ,6)*y.^2;

        %Obtain local error estimate with quadrature
        interior= @(x,y) (u_cell(x,y) - exact_sol(x,y)).^2;
        local_error =int_over_cell(interior,cell,nod_in_e2d,tri,xcoord,ycoord); %Obtain integral via  Legendre-Gauss quadrature
        global_L2_Error = global_L2_Error + local_error;
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


function cell_avg = build_cell_avg_vec(inner_v_vec,f, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area)
    cell_avg = NaN(1,n2d);
    for ii = 1:length(inner_v_vec)
        i = inner_v_vec(ii);

        cell_int = int_over_cell(f,i,nod_in_e2d,tri,xcoord,ycoord);
        cell_avg(i) = cell_int/area(i);
    end
end


function ref_avg = build_ref_array(C,area,n_steps, ind_vec,nod_in_e2d,tri,xcoord,ycoord,n2d)
    ref_avg = NaN(1,n2d);
    for ii = 1:n_steps
        i = ind_vec(ii);
        f = @(x,y) C(ii,1) + C(ii ,2)*x + C(ii ,3)*y+ C(ii ,4)*x.*y + C(ii ,5)*x.^2 + C(ii ,6)*y.^2;
        cell_int = int_over_cell(f,i,nod_in_e2d,tri,xcoord,ycoord);
        ref_avg(i) = cell_int/area(i);
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
    
    cell_int = cell_int + tri_int(lil_t_1, f) +  tri_int(lil_t_2, f);
    end
end

function result = tri_int(tri_coord, f)

%Updated to be symmetric, should work polynomials degree 2


%tri_coord : 3 x 2 containing triangle coordinates to integrate over

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


function mass_m = generate_mass_matrix(nghbr_n2d, n2d, e2d,elem_area,tri,ng_num)
  mass_m=zeros(size(nghbr_n2d));
  aux=zeros([1,n2d]);
  for j=1:e2d,
      enodes=tri(:,j);
      for n=1:3,
      row=enodes(n);
      cc=ng_num(row);
      aux(nghbr_n2d(1:cc,row))=1:cc;
         for m=1:3,
            col=enodes(m);
	    pos=aux(col);
	    mass_m(pos,row)=mass_m(pos,row)+0.25*elem_area(j)/3.0;
	    if m==n,  
            mass_m(pos,row)=mass_m(pos,row)+0.25*elem_area(j)/3.0;
	    end;
	  end;   
       end;
  end;     

    % correct mass matrix so that it agrees with nghbr_n2d 
  mass_m=zeros(size(nghbr_n2d));
  aux=zeros([1,n2d]);
  for j=1:e2d,
      enodes=tri(:,j);
      for n=1:3,
      row=enodes(n);
      cc=ng_num(row);
      aux(nghbr_n2d(1:cc,row))=1:cc;
         for m=1:3,
            col=enodes(m);
	    pos=aux(col);
	    mass_m(pos,row)=mass_m(pos,row)+0.25*elem_area(j)/3.0;
	    if m==n,  
            mass_m(pos,row)=mass_m(pos,row)+0.25*elem_area(j)/3.0;
	    end;
	  end;   
      end;
  end;
 end


 function [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d)

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




