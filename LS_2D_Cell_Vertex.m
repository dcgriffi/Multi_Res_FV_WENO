%LS_2D_Cell_Vertex.m
%   Constrained least squares reconstruction of a quadratic polynomial
%   (k = 6 coefficients: 1, x, y, xy, x^2, y^2) from cell averages on the
%   median-dual control volumes of a regular triangular mesh. The stencil
%   is the target control volume and its six neighbors (m = 7). The
%   polynomial matches the target cell average exactly and the remaining
%   cells in the least squares sense.
%
%   test = "example"     : coarse mesh of [0,1]^2, f = sin(2*pi*(x+y)). For
%                          one target vertex, prints the largest difference
%                          of the cell averages of p and f on the stencil and
%                          saves a figure of the differences on the control
%                          volumes to Results/LS_2D/ls_hex_example.pdf.
%   test = "convergence" : L2 error of the reconstruction of sin(2*pi*(x+y))
%                          on meshes of [0,1]^2 with spacing h = 1/n.
%                          Expected order: 3.
%
%   The least squares system is built in local coordinates centered at the
%   vertex and scaled by sqrt(area), then mapped back to 1, x, y, ...
%   Mesh routines (gen_reg_tri_mesh.m, get_mesh_sten_info) are modified from
%   code by Sergey Danilov.

clear; close all;

test = "example";   % "example" or "convergence"

k = 6; m = 7;            % quadratic polynomial; vertex and its six neighbors
f_exact = @(x,y) sin(2*pi*(x+y));
xL = 0; xR = 1;

if test == "example"
    n_array = 6;
else
    n_array = [10,20,40,80,160];
end
h_array = zeros(length(n_array),1); L2_Error = zeros(length(n_array),1);

for kk = 1:length(n_array)
    n = n_array(kk); h = (xR - xL)/n; h_array(kk) = h;

    %Regular triangular mesh of [0,1]^2 padded with two layers of vertices
    x = (xL:h:xR); y = x;
    x = [x(1) - 2*h, x(1) - h, x, x(end) + h, x(end) + 2*h];
    y = [y(1) - 2*h, y(1) - h, y, y(end) + h, y(end) + 2*h];
    N = n + 5;                                   % vertices per direction
    [tri, nodnum, xcoord, ycoord, e2d, n2d] = gen_reg_tri_mesh(x,y);
    [nod_in_e2d, nghbr_n2d, elem_area, area, n_num, ng_num, edges, ed2d] = get_mesh_sten_info(tri, xcoord, ycoord, e2d, n2d);

    %Vertices with complete control volumes, and those whose whole stencil
    %has complete control volumes (where we reconstruct)
    all_v = reshape(1:N*N, N, N);
    inner_v_vec = reshape(all_v(2:end-1, 2:end-1), [], 1);
    ind_vec = reshape(all_v(3:end-2, 3:end-2), [], 1);

    cell_avg = build_cell_avg_vec(inner_v_vec, f_exact, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);

    if test == "example"
        i = all_v(ceil(N/2), ceil(N/2));        % target vertex near the center
        stencil = nghbr_n2d(:,i);
        W = hexagon_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area);
        C = W*[reshape(cell_avg(stencil(1:m)), [], 1); -cell_avg(i)];
        p = @(xx,yy) C(1) + C(2)*xx + C(3)*yy + C(4)*xx.*yy + C(5)*xx.^2 + C(6)*yy.^2;
        p_avg = build_cell_avg_vec(inner_v_vec, p, tri, nod_in_e2d, n2d, e2d, xcoord, ycoord, area);
        d = p_avg - cell_avg;
        fprintf('Largest |p_avg - f_avg| on the stencil: %.3e\n', max(abs(d(stencil(1:m)))));
        pth = fullfile(pwd, 'Results', 'LS_2D');
        if ~exist(pth, 'dir'), mkdir(pth); end
        plot_cells = reshape(all_v(4:end-3, 4:end-3), [], 1);   % drop one outer layer for the figure
        plot_hex_grid(d, plot_cells, i, stencil(1:m), nod_in_e2d, tri, xcoord, ycoord, ...
            fullfile(pth, 'ls_hex_example.pdf'));
        return
    end

    %Reconstruct in every vertex of ind_vec and integrate the squared error
    err2 = 0;
    for ii = 1:length(ind_vec)
        i = ind_vec(ii);
        stencil = nghbr_n2d(:,i);
        W = hexagon_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area);
        C = W*[reshape(cell_avg(stencil(1:m)), [], 1); -cell_avg(i)];
        p = @(xx,yy) C(1) + C(2)*xx + C(3)*yy + C(4)*xx.*yy + C(5)*xx.^2 + C(6)*yy.^2;
        err2 = err2 + int_over_cell(@(xx,yy) (p(xx,yy) - f_exact(xx,yy)).^2, i, nod_in_e2d, tri, xcoord, ycoord);
    end
    L2_Error(kk) = sqrt(err2);
end

order = [NaN; log(L2_Error(1:end-1)./L2_Error(2:end))./log(h_array(1:end-1)./h_array(2:end))];
disp(table(n_array(:), L2_Error, order, 'VariableNames', {'n','L2 error','L2 order'}))


function plot_hex_grid(d, cells, i, stencil, nod_in_e2d, tri, xcoord, ycoord, fname)
    %Median-dual control volumes of the vertices in cells, with the value of
    %d in each. Stencil cells are shaded light gray and their outer boundary
    %is outlined in red; the target vertex i is shaded darker. No title;
    %details go in the LaTeX caption.
    fig = figure('Units','inches','Position',[1 1 5.5 5.0],'Color','w'); hold on
    E = zeros(0,4);                               % edges of the stencil cells
    for c = cells(:)'
        P = dual_cell_polygon(c, nod_in_e2d, tri, xcoord, ycoord);
        if c == i
            col = [0.75 0.75 0.75];
        elseif any(stencil == c)
            col = [0.92 0.92 0.92];
        else
            col = [1 1 1];
        end
        patch(P(:,1), P(:,2), col, 'EdgeColor', [0.5 0.5 0.5]);
        text(xcoord(c), ycoord(c), sprintf('%.1e', d(c)), 'HorizontalAlignment','center', 'FontSize', 9);
        if any(stencil == c), E = [E; P, circshift(P, -1)]; end %#ok<AGROW>
    end
    %Outer boundary of the stencil: edges that belong to only one stencil cell
    key = round(sort([E(:,1)+1i*E(:,2), E(:,3)+1i*E(:,4)], 2, 'ComparisonMethod', 'real')*1e8);
    [~, ~, g] = unique([real(key), imag(key)], 'rows');
    once = accumarray(g, 1) == 1;
    for q = find(once(g))'
        plot([E(q,1) E(q,3)], [E(q,2) E(q,4)], 'r-', 'LineWidth', 2);
    end
    axis equal; axis off;
    exportgraphics(fig, fname, 'ContentType', 'vector');
end

function P = dual_cell_polygon(v, nod_in_e2d, tri, xcoord, ycoord)
    %Vertices of the median-dual control volume of vertex v: the centroids
    %of the surrounding triangles and the midpoints of the edges at v,
    %ordered by angle around v.
    pts = [];
    for t = nod_in_e2d(:,v)'
        if t == 0, continue; end
        vt = tri(:,t);
        pts = [pts; mean(xcoord(vt)), mean(ycoord(vt))]; %#ok<AGROW>
        for w = vt(vt ~= v)'
            pts = [pts; (xcoord(v) + xcoord(w))/2, (ycoord(v) + ycoord(w))/2]; %#ok<AGROW>
        end
    end
    pts = uniquetol(pts, 1e-12, 'ByRows', true);
    [~, o] = sort(atan2(pts(:,2) - ycoord(v), pts(:,1) - xcoord(v)));
    P = pts(o,:);
end

function W = hexagon_condensed_constrained_weight_matrix(k,m,i,nghbr_n2d, nod_in_e2d,tri,xcoord,ycoord,area)

    
        %Center the least squares on vertex i and scale by sqrt(area). Building
        %the system in absolute coordinates makes A'A badly conditioned.
        xc = xcoord(i); yc = ycoord(i); s = sqrt(area(i));
        f1 = @(x,y) 1;    f2 = @(x,y) (x-xc)/s;       f3 = @(x,y) (y-yc)/s; 
        f4 = @(x,y) ((x-xc)/s).*((y-yc)/s); f5 = @(x,y) ((x-xc)/s).^2; f6 = @(x,y) ((y-yc)/s).^2;
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

        %Map coefficients of 1, z, w, zw, z^2, w^2 (z = (x-xc)/s, w = (y-yc)/s)
        %back to coefficients of 1, x, y, xy, x^2, y^2, so the rest of the
        %code is unchanged.
        W = scaled_to_absolute_2d(xc,yc,s)*W;



end


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
