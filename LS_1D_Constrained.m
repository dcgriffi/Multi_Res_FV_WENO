%LS_1D_Constrained.m
%   Constrained least squares reconstruction of a quadratic polynomial from
%   cell averages in one dimension (k = 3 coefficients, m = 5 cells, r = 2
%   cells to the left of the target cell). The polynomial matches the cell
%   average of the target cell exactly and the remaining stencil cells in
%   the least squares sense. This is the building block of the
%   multi-resolution WENO reconstruction in Multires_1D.m.
%
%   test = "example"     : N = 10 cells on [0,1], f = sin(2*pi*x), target
%                          cell i = 3. Prints the cell averages of f and of
%                          the polynomial on the stencil and saves a plot to
%                          Results/LS_1D/ls_constrained_example.pdf.
%   test = "convergence" : L2 error of the reconstruction of sin(2*pi*x) on
%                          a uniform periodic grid of [0,1], N = 20,...,640.
%                          Expected order: 3.
%
%   The constrained problem is solved with a Lagrange multiplier. Its
%   solution can be written as C = W*[d; -u_i], where d holds the stencil
%   cell averages and u_i is the target cell average, so W is precomputed
%   once per cell. The system is built in local coordinates
%   z = (x - x_c)/h_c centered at the target cell and scaled by its width.
%   In the original coordinates, (A'*A)\L' loses accuracy on fine meshes
%   (the reconstruction error stagnates near N = 640 on [0,1]). The
%   coefficients are mapped back to the monomial basis 1, x, x^2.
%
%   Requires lgwt.m.

clear; close all;

test = "convergence";   % "example" or "convergence"

k = 3; m = 5; r = 2;    % quadratic polynomial, five-cell centered stencil
f = @(x) sin(2*pi*x);
xL = 0; xR = 1;

pth = fullfile(pwd, 'Results', 'LS_1D');    % output folder (ignored by git)
if ~exist(pth, 'dir'), mkdir(pth); end

if test == "example"
    N = 10; h = (xR-xL)/N; x = (xL:h:xR)';
    i = 3;                                      % target cell [0.2, 0.3]

    fbar = cell_averages(f, x);
    d = [fbar(i-2:i+2); -fbar(i)];
    C = get_constrained_weight_matrix(k,m,r,i,x)*d;
    p = @(z) C(1) + C(2)*z + C(3)*z.^2;
    pbar = cell_averages(p, x);

    T = table((i-2:i+2)', x(i-2:i+2), x(i-1:i+3), fbar(i-2:i+2), pbar(i-2:i+2), ...
        pbar(i-2:i+2) - fbar(i-2:i+2), 'VariableNames', ...
        {'j','x_left','x_right','f_avg','p_avg','p_avg - f_avg'});
    disp(T)

    %Publication-style plot (no title; details go in the figure caption)
    fig = figure('Units','inches','Position',[1 1 5.5 4.0],'Color','w'); hold on
    patch([x(i-2) x(i+3) x(i+3) x(i-2)], [-1.2 -1.2 1.2 1.2], [0.92 0.92 0.92], 'EdgeColor','none');
    patch([x(i) x(i+1) x(i+1) x(i)], [-1.2 -1.2 1.2 1.2], [0.80 0.80 0.80], 'EdgeColor','none');
    for j = 1:N
        hc = plot([x(j) x(j+1)], [fbar(j) fbar(j)], 'k-', 'LineWidth', 2);
    end
    xs = linspace(x(i-2), x(i+3), 200);
    hp = plot(xs, p(xs), 'r-', 'LineWidth', 1.5);
    xf = linspace(xL, xR, 400);
    he = plot(xf, f(xf), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    xlim([xL xR]); ylim([-1.2 1.2]); xlabel('x'); ylabel('u');
    legend([hc hp he], {'Cell averages of f','Constrained least squares p', ...
        'f(x) = sin(2\pi x)'}, 'Location','southwest','Box','off');
    set(gca,'FontSize',13,'LineWidth',1,'TickDir','out','Box','on','Layer','top');
    exportgraphics(fig, fullfile(pth, 'ls_constrained_example.pdf'), 'ContentType', 'vector');

elseif test == "convergence"
    n_array = [20,40,80,160,320,640];
    dx = zeros(length(n_array),1); L2_Error = zeros(length(n_array),1);

    for kk = 1:length(n_array)
        N = n_array(kk);
        dx(kk) = (xR-xL)/N; h = dx(kk);
        x = (xL:h:xR)';
        u_avg = cell_averages(f, x);

        %Coefficients in every cell (periodic data)
        up2 = circshift(u_avg,-2); up1 = circshift(u_avg,-1);
        um1 = circshift(u_avg, 1); um2 = circshift(u_avg, 2);
        C = zeros(N,k);
        for i = 1:N
            W = get_constrained_weight_matrix(k,m,r,i,x);
            C(i,:) = W*[um2(i); um1(i); u_avg(i); up1(i); up2(i); -u_avg(i)];
        end

        L2_Error(kk) = sqrt(global_L2_error(C,f,N,x));
    end
    order_table(n_array, dx, L2_Error);
end


function W = get_constrained_weight_matrix(k,m,r,i,x_array)
%Weight matrix W with C = W*[d; -u_i] for the constrained least squares
%problem: minimize |A*C - d| subject to L*C = u_i. Row j of A holds the cell
%averages of 1, x, x^2 over the j-th stencil cell and L holds them over the
%target cell. Built in scaled local coordinates, then mapped back.
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
L = zeros(1,k);
for l = 1:k
    L(l) = (z(mod_n(i + 1))^l - z(mod_n(i))^l)/(l * dz(mod_n(i)));
end

%Schur complement solution of the saddle point system
[Q,R] = qr(A,0);
W = R\(Q');                                  % (A'*A)\A'
P = R\(R'\L');                               % (A'*A)\L'
mu = -L*P;
W = [W + (P*L*W)/mu, P/mu];

%Map coefficients of 1, z, z^2 back to coefficients of 1, x, x^2
T = [1, -xc/hc, xc^2/hc^2;
     0,  1/hc, -2*xc/hc^2;
     0,  0,     1/hc^2];
W = T*W;
end

function avg = cell_averages(fun, x)
%Cell averages of fun over the cells [x(j), x(j+1)] (5-point Gauss)
N = length(x) - 1; avg = zeros(N,1);
for j = 1:N
    [vals,weights] = lgwt(5,x(j),x(j+1));
    avg(j) = sum(fun(vals).*weights)/(x(j+1) - x(j));
end
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
