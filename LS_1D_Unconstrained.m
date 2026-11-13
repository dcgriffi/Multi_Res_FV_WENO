%LS_1D_Unconstrained.m
%   Unconstrained least squares reconstruction of a quadratic polynomial
%   from cell averages in one dimension (k = 3 coefficients, m = 3 cells,
%   r = 1 cell to the left of the target cell).
%
%   test = "check"       : non-uniform grid with arbitrary cell averages.
%                          With m = k the least squares problem is square,
%                          so the polynomial matches all three cell averages
%                          exactly (printed differences are at round-off).
%   test = "convergence" : L2 error of the reconstruction of sin(2*pi*x) on
%                          a uniform periodic grid of [0,1], N = 20,...,640.
%                          Expected order: 3.
%
%   The least squares system is built in local coordinates
%   z = (x - x_c)/h_c centered at the target cell and scaled by its width,
%   which keeps the system well conditioned as the mesh is refined. The
%   coefficients are then mapped back to the monomial basis 1, x, x^2.
%
%   Requires lgwt.m.

clear; close all;

test = "convergence";   % "check" or "convergence"

k = 3; m = 3; r = 1;    % quadratic polynomial, three-cell centered stencil

if test == "check"
    %Non-uniform grid and arbitrary cell averages
    x_array = [0; 0.5; 0.9; 1.4; 2.3; 3.2; 4];
    v_array = [-1; 3; 2; 4; 0; 1];
    dx = x_array(2:end) - x_array(1:end-1);
    n = length(x_array);
    i = 3;                                      % target cell

    W = get_weight_matrix(k,m,r,i,x_array);
    d = [v_array(i-1); v_array(i); v_array(i+1)];
    C = W*d;
    q = @(x) C(1) + C(2)*x + C(3)*x.^2;

    %Cell averages of q over the stencil cells
    q_avg = zeros(n-1,1);
    for j = 1:n-1
        [vals,weights] = lgwt(5,x_array(j),x_array(j+1));
        q_avg(j) = sum(q(vals).*weights)/dx(j);
    end
    disp('Difference of cell averages (q - data) on the stencil cells i-1, i, i+1:')
    disp((q_avg(i-1:i+1) - v_array(i-1:i+1))')

elseif test == "convergence"
    xL = 0; xR = 1;
    u_exact = @(x) sin(2*pi*x);
    n_array = [20,40,80,160,320,640];
    dx = zeros(length(n_array),1); L2_Error = zeros(length(n_array),1);

    for kk = 1:length(n_array)
        N = n_array(kk);
        dx(kk) = (xR-xL)/N; h = dx(kk);
        x = (xL:h:xR)';

        %Exact cell averages
        u_avg = zeros(N,1);
        for j = 1:N
            u_avg(j) = integral(u_exact, x(j), x(j+1), 'ArrayValued', true)/h;
        end

        %Coefficients in every cell (periodic data)
        up1 = circshift(u_avg,-1); um1 = circshift(u_avg,1);
        C = zeros(N,k);
        for i = 1:N
            W = get_weight_matrix(k,m,r,i,x);
            C(i,:) = W*[um1(i); u_avg(i); up1(i)];
        end

        L2_Error(kk) = sqrt(global_L2_error(C,u_exact,N,x));
    end
    order_table(n_array, dx, L2_Error);
end


function W = get_weight_matrix(k,m,r,i,x_array)
%Weight matrix W with C = W*d for the unconstrained least squares problem
%A*C = d, where row j of A holds the cell averages of 1, x, x^2 over the
%j-th stencil cell. Built in scaled local coordinates, then mapped back.
x_last_ghost = x_array(end) + x_array(2:k + 2);
x_first_ghost = x_array(1) - (x_array(end) - x_array(end - (k - 1):end-1));
x_ghost = [x_first_ghost; x_array; x_last_ghost];
dx_ghost = x_ghost(2:end) - x_ghost(1:end-1);
mod_n = @(i) i + k - 1;

%Local coordinates centered at the target cell and scaled by its width
xc = (x_ghost(mod_n(i)) + x_ghost(mod_n(i+1)))/2;
hc = dx_ghost(mod_n(i));
z = (x_ghost - xc)/hc; dz = dx_ghost/hc;

A = zeros(m,k);
for l = 1:k
    for j = 1:m
        A(j,l) = (z(mod_n(i - r + j))^l - z(mod_n(i - r + j - 1))^l)/(l * dz(mod_n(i - r + j - 1)));
    end
end

[Q,R] = qr(A,0);
W = R\(Q');                                  % (A'*A)\A' without forming A'*A

%Map coefficients of 1, z, z^2 back to coefficients of 1, x, x^2
T = [1, -xc/hc, xc^2/hc^2;
     0,  1/hc, -2*xc/hc^2;
     0,  0,     1/hc^2];
W = T*W;
end

function err2 = global_L2_error(C,u_exact,N,x)
%Squared L2 error of the piecewise quadratic reconstruction (5-point Gauss)
err2 = 0;
for i = 1:N
    [vals,weights] = lgwt(5,x(i),x(i+1));
    q = C(i,1) + C(i,2)*vals + C(i,3)*vals.^2;
    err2 = err2 + sum((q - u_exact(vals)).^2.*weights);
end
end

function order_table(n_array, dx, err)
order = [NaN; log(err(1:end-1)./err(2:end))./log(dx(1:end-1)./dx(2:end))];
T = table(n_array(:), err(:), order, 'VariableNames', {'N','L2 error','L2 order'});
disp(T)
end
