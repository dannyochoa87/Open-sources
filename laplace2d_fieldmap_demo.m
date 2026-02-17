% laplace2d_fieldmap_demo.m
% ------------------------------------------------------------
% Universidad de Cuenca
%
% Purpose (paper-aligned workflow):
%   This script reproduces the numerical pipeline described in the article:
%   (1) Define geometry and Dirichlet boundary conditions:
%       - Outer grounded contour (V = V0) and/or conductor equipotentials (V = Vm)
%       - Fixed-node mask chi_{i,j} that marks Dirichlet nodes
%       (paper: Eq. \eqref{eq:dirichlet_bc}--\eqref{eq:conductor_dirichlet}, Eq. \eqref{eq:mask})
%   (2) Discretize the 2D Laplace problem on a rectangular grid:
%       - V_{i,j} approximates V(x_i, y_j)
%       (paper: Eq. \eqref{eq:laplace}, Eq. \eqref{eq:grid_nodes})
%   (3) Solve the discrete Laplace equation by fixed-point iterations on free nodes:
%       - Update rule from the finite-difference Laplacian
%       (paper: Eq. \eqref{eq:update_general} or Eq. \eqref{eq:update_uniform})
%       - Stopping rule based on max voltage change between iterations
%       (paper: Eq. \eqref{eq:stop_delta})
%   (4) Post-process the electric field:
%       - E = -grad(V) and finite-difference estimates for Ex, Ey
%       (paper: Eq. \eqref{eq:E_def}, Eq. \eqref{eq:Ex_disc}--\eqref{eq:Ey_disc})
%   (5) Plot V maps, equipotentials, and E visualizations.
%
% Open-code use and modifications are permitted, provided that the
% scientific article documenting this tool is properly cited.
% ------------------------------------------------------------

clear; clc; close all;

fprintf('\n2D Laplace Electrostatic Mapper (FDM + fixed-point)\n');
fprintf('--------------------------------------------------\n');
fprintf('Cases:\n');
fprintf('  1) Tripolar cable (grounded outer screen)\n');
fprintf('  2) Trench installation (grounded rectangular boundary)\n');
% fprintf('  3) Rectangular box (Dirichlet boundary only; no conductors)\n');
% fprintf('  4) Coaxial-like layout (3 conductors + grounded outer circle)\n\n');

caseID = askInt('Select a case [1-2]', 1, 4, 1);

% --------------------------- Solver parameters ---------------------------
% eps  : voltage-change stopping tolerance (paper: Eq. \eqref{eq:stop_delta})
% kmax : hard cap on iteration count (practical safeguard)
P = struct();
P.eps  = askNum('Tolerance eps (V) for ||ΔV||_inf', 1e-2);
P.kmax = askInt('Max iterations kmax', 10, 200000, 6000);

% Default excitation snapshot (used as a simple reproducible example)
Va_def = 100; Vb_def = 50; Vc_def = -10;

% --------------------------- Case-specific inputs ---------------------------
% These blocks define:
%   - the computational window and geometry parameters
%   - the Dirichlet values on conductors and grounded contours
% The Dirichlet information is enforced via chi_{i,j} (paper: Eq. \eqref{eq:mask})
switch caseID
    case 1
        % Tripolar cable: circular region inside an outer grounded screen
        % Outer boundary (screen): V = V0 (paper: Eq. \eqref{eq:dirichlet_bc})
        % Conductors: V = VA, VB, VC on internal equipotential regions (paper: Eq. \eqref{eq:conductor_dirichlet})
        P.r1    = askNum('Outer screen radius r1 (m)', 0.25);
        P.rcond = askNum('Conductor radius r_cond (m)', 0.05);
        P.rc    = askNum('Radius of conductor-center circle r_c (m)', 0.10);
        P.V0    = askNum('Ground/screen potential V0 (V)', 0);
        P.VA    = askNum('Phase A potential VA (V)', Va_def);
        P.VB    = askNum('Phase B potential VB (V)', Vb_def);
        P.VC    = askNum('Phase C potential VC (V)', Vc_def);

        defaultMesh = [10 20 50 100];

    case 2
        % Trench: rectangular boundary is grounded (Dirichlet: V = V0 on all edges)
        % Conductors are circular internal equipotentials
        P.A     = askNum('Trench width A (m)', 0.50);
        P.B     = askNum('Trench height B (m)', 0.50);
        P.C     = askNum('Burial depth from top C (m)', 0.25);    % conductor center y = B - C
        P.D     = askNum('Half horizontal spacing D (m)', 0.125); % centers at x = A/2 +- D and A/2
        P.rcond = askNum('Conductor radius r_cond (m)', 0.05);
        P.V0    = askNum('Ground potential V0 (V)', 0);
        P.VA    = askNum('Phase A potential VA (V)', Va_def);
        P.VB    = askNum('Phase B potential VB (V)', Vb_def);
        P.VC    = askNum('Phase C potential VC (V)', Vc_def);

        defaultMesh = [16 32 64 128];

    case 3
        % Rectangular box: Dirichlet boundaries only (no conductors)
        % Useful as a baseline Laplace example.
        P.Lx = askNum('Domain width Lx', 1.0);
        P.Ly = askNum('Domain height Ly', 1.0);

        useConst = askYesNo('Constant boundary on all sides?', true);
        P.boxConst = useConst;

        if useConst
            P.V0 = askNum('Boundary potential V0 (V)', 0);
        else
            P.Vtop    = askNum('Top boundary Vtop (V)', 0);
            P.Vbottom = askNum('Bottom boundary Vbottom (V)', 0);
            P.Vleft   = askNum('Left boundary Vleft (V)', 0);
            P.Vright  = askNum('Right boundary Vright (V)', 0);
        end

        defaultMesh = [20 50 100];

    case 4
        % Coaxial-like: circular grounded contour plus three internal conductors
        P.Rout  = askNum('Outer grounded radius Rout (m)', 0.50);
        P.rcent = askNum('Radius of conductor-center circle (m)', 0.25);
        P.rcond = askNum('Conductor radius r_cond (m)', 0.08);
        P.V0    = askNum('Outer boundary potential V0 (V)', 0);
        P.VA    = askNum('Conductor A potential VA (V)', Va_def);
        P.VB    = askNum('Conductor B potential VB (V)', Vb_def);
        P.VC    = askNum('Conductor C potential VC (V)', Vc_def);

        defaultMesh = [33 65 129];
end

meshList = askMeshList('Mesh sizes N (space-separated). Each N gives N×N nodes.', defaultMesh);

fprintf('\nRunning meshes: ');
fprintf('%d ', meshList);
fprintf('\n\n');

% --------------------------- Sweep over mesh resolutions ---------------------------
% For each N:
%   1) Build grid (paper: Eq. \eqref{eq:grid_nodes})
%   2) Build chi_{i,j} mask and set Dirichlet values on chi=1 nodes (paper: Eq. \eqref{eq:mask})
%   3) Iterate fixed-point update on chi=0 nodes (paper: Eq. \eqref{eq:update_general})
%   4) Stop when ||ΔV||_inf < eps (paper: Eq. \eqref{eq:stop_delta})
%   5) Compute E = -grad(V) (paper: Eq. \eqref{eq:E_def})
Res = struct([]);
for m = 1:numel(meshList)
    N = meshList(m);
    t0 = tic;

    % (1) Grid construction
    % x_i, y_j spacing hx, hy define the FDM mesh (paper: Eq. \eqref{eq:grid_nodes})
    [x, y, X, Y, hx, hy] = buildGrid(caseID, P, N);

    % (2) Dirichlet sets and fixed-node mask chi_{i,j}
    % chi = true  -> fixed node (boundary or conductor)
    % chi = false -> free node (unknown updated by iterations)
    % (paper: Eq. \eqref{eq:mask})
    [V0, chi, plotNaNMask] = buildMaskAndDirichlet(caseID, P, X, Y);

    % (3) Fixed-point iterations on free nodes (chi=0)
    % Discrete Laplacian update (paper: Eq. \eqref{eq:update_general})
    % Stopping rule (paper: Eq. \eqref{eq:stop_delta})
    [Vsol, iters, dmax, hist] = laplaceFixedPoint(V0, chi, hx, hy, P.eps, P.kmax);

    % (4) Electric field post-processing
    % E(x,y) = -∇V (paper: Eq. \eqref{eq:E_def})
    % Ex, Ey via finite differences (paper: Eq. \eqref{eq:Ex_disc}--\eqref{eq:Ey_disc})
    [Ex, Ey] = electricField(Vsol, hx, hy);
    Eabs = hypot(Ex, Ey); % |E| (paper: Eq. \eqref{eq:E_mag})

    % Store results
    Res(m).N = N;
    Res(m).x = x; Res(m).y = y;
    Res(m).X = X; Res(m).Y = Y;
    Res(m).V = Vsol;
    Res(m).Ex = Ex; Res(m).Ey = Ey; Res(m).Eabs = Eabs;
    Res(m).iters = iters;
    Res(m).dmax  = dmax;
    Res(m).time_s = toc(t0);
    Res(m).hist = hist;
    Res(m).plotNaNMask = plotNaNMask;
end

% --------------------------- Summary table (command window) ---------------------------
fprintf('Summary (convergence by ||ΔV||_inf):\n');
fprintf('  N×N\titers\tfinal ||ΔV||_inf (V)\ttime (s)\n');
for m = 1:numel(Res)
    fprintf('  %d\t%d\t%.4e\t\t\t%.3f\n', Res(m).N, Res(m).iters, Res(m).dmax, Res(m).time_s);
end
fprintf('\n');

% --------------------------- Convergence traces ---------------------------
% This figure visualizes the evolution of ||ΔV||_inf vs iteration k
% (paper: Eq. \eqref{eq:stop_delta})
if numel(Res) > 1
    figure('Name','Convergence traces (||ΔV||_inf)','Color','w');
    hold on; grid on;
    for m = 1:numel(Res)
        semilogy(Res(m).hist, 'LineWidth', 1.2);
    end
    xlabel('Iteration k');
    ylabel('||\Delta V||_\infty (V)');
    legend(arrayfun(@(r)sprintf('%dx%d', r.N, r.N), Res, 'UniformOutput', false), 'Location','northeast');
    title('Fixed-point convergence');
    hold off;
end

% --------------------------- Plots for finest mesh ---------------------------
% The last entry is assumed to be the finest mesh (largest N).
% Plots:
%   - V heat map (imagesc)
%   - equipotential contours (contourf)
%   - E vectors (quiver) over faint equipotentials
%   - |E| map (imagesc)
R = Res(end);
Vplot = R.V;

% For circular domains, hide the outside-of-domain region in plots
if ~isempty(R.plotNaNMask)
    Vplot(R.plotNaNMask) = NaN;
end

% (A) Potential heat map
figure('Name','Potential heatmap','Color','w');
imagesc(R.x, R.y, Vplot);
set(gca,'YDir','normal'); axis equal tight;
cb = colorbar; cb.Label.String = 'V (V)';
colormap(getCmap());
xlabel('x'); ylabel('y');
title(sprintf('Potential heatmap (%dx%d)', R.N, R.N));

% (B) Equipotential contours
figure('Name','Equipotentials','Color','w');
contourf(R.X, R.Y, Vplot, 20, 'LineWidth', 0.5);
set(gca,'YDir','normal'); axis equal tight;
cb = colorbar; cb.Label.String = 'V (V)';
colormap(getCmap());
xlabel('x'); ylabel('y');
title(sprintf('Equipotential contours (%dx%d)', R.N, R.N));

% (C) Electric field vectors over faint equipotentials
figure('Name','Electric field vectors','Color','w');
hold on; grid on;
contour(R.X, R.Y, Vplot, 12, 'LineColor', [0.75 0.75 0.75], 'LineWidth', 0.8);
step = max(1, round(R.N/20)); % downsample vectors for readability
quiver(R.X(1:step:end,1:step:end), R.Y(1:step:end,1:step:end), ...
       R.Ex(1:step:end,1:step:end), R.Ey(1:step:end,1:step:end), ...
       1.5, 'LineWidth', 0.9);
set(gca,'YDir','normal'); axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('E-field vectors (%dx%d)', R.N, R.N));
hold off;

% (D) |E| magnitude map
figure('Name','|E| magnitude','Color','w');
Eplot = R.Eabs;
if ~isempty(R.plotNaNMask)
    Eplot(R.plotNaNMask) = NaN;
end
imagesc(R.x, R.y, Eplot);
set(gca,'YDir','normal'); axis equal tight;
cb = colorbar; cb.Label.String = '|E| (V/m)';
colormap(getCmap());
xlabel('x'); ylabel('y');
title(sprintf('|E| magnitude (%dx%d)', R.N, R.N));

% Optional: print an iteration-count table for trench sweeps
if caseID == 2 && numel(Res) > 1
    fprintf('Iteration counts vs mesh (trench case):\n');
    fprintf('  N×N\titers\n');
    for m = 1:numel(Res)
        fprintf('  %d\t%d\n', Res(m).N, Res(m).iters);
    end
    fprintf('\n');
end

% ======================================================================
% Local functions
% ======================================================================

function val = askNum(prompt, defaultVal)
% Reads a numeric value from the command window with a default fallback.
    s = input(sprintf('%s [default = %g]: ', prompt, defaultVal), 's');
    if isempty(strtrim(s))
        val = defaultVal;
    else
        val = str2double(s);
        if isnan(val)
            warning('Invalid input. Using default.');
            val = defaultVal;
        end
    end
end

function val = askInt(prompt, minVal, maxVal, defaultVal)
% Reads an integer value from the command window with bounds and default.
    s = input(sprintf('%s [default = %d]: ', prompt, defaultVal), 's');
    if isempty(strtrim(s))
        val = defaultVal;
        return;
    end
    val = round(str2double(s));
    if isnan(val) || val < minVal || val > maxVal
        warning('Invalid input. Using default.');
        val = defaultVal;
    end
end

function tf = askYesNo(prompt, defaultTF)
% Reads a yes/no answer from the command window.
    def = 'y';
    if ~defaultTF, def = 'n'; end
    s = lower(strtrim(input(sprintf('%s [y/n, default = %s]: ', prompt, def), 's')));
    if isempty(s), tf = defaultTF; return; end
    tf = startsWith(s,'y');
end

function meshList = askMeshList(prompt, defaultList)
% Reads a list of mesh sizes N; each N implies an N×N grid.
    defStr = sprintf('%d ', defaultList);
    s = input(sprintf('%s\n  Enter N list [default = %s]: ', prompt, strtrim(defStr)), 's');
    if isempty(strtrim(s))
        meshList = defaultList(:).';
        return;
    end
    nums = str2num(s); %#ok<ST2NM>
    if isempty(nums) || any(nums < 5)
        warning('Invalid mesh list. Using default.');
        meshList = defaultList(:).';
    else
        meshList = unique(round(nums(:).'), 'stable');
    end
end

function [x, y, X, Y, hx, hy] = buildGrid(caseID, P, N)
% Builds the computational grid (paper: Eq. \eqref{eq:grid_nodes}).
% Outputs:
%   x, y : 1D node vectors
%   X, Y : meshgrid arrays
%   hx, hy : uniform spacings
    switch caseID
        case 1
            % Square window enclosing the cable circle
            x = linspace(-P.r1, P.r1, N);
            y = linspace(-P.r1, P.r1, N);
        case 2
            x = linspace(0, P.A, N);
            y = linspace(0, P.B, N);
        case 3
            x = linspace(0, P.Lx, N);
            y = linspace(0, P.Ly, N);
        case 4
            x = linspace(-P.Rout, P.Rout, N);
            y = linspace(-P.Rout, P.Rout, N);
    end
    [X, Y] = meshgrid(x, y);
    hx = x(2) - x(1);
    hy = y(2) - y(1);
end

function [V, chi, plotNaNMask] = buildMaskAndDirichlet(caseID, P, X, Y)
% Creates:
%   V   : initial voltage matrix with Dirichlet values placed at fixed nodes
%   chi : fixed-node mask (true for fixed nodes, false for free nodes)
%         (paper: Eq. \eqref{eq:mask})
%   plotNaNMask : optional mask to hide outside-of-domain regions in plots
%
% Dirichlet enforcement matches:
%   - external boundary: V = \bar V(x,y) on ∂Ω (paper: Eq. \eqref{eq:dirichlet_bc})
%   - conductors: V = \bar V_m on Γ_m (paper: Eq. \eqref{eq:conductor_dirichlet})
    V = zeros(size(X));
    chi = false(size(X));
    plotNaNMask = [];

    switch caseID
        case 1
            % Tripolar cable with grounded outer screen
            r = hypot(X, Y);
            outer = (r >= P.r1); % treat outside as fixed (screen reference)

            % Conductor centers placed on a circle (trefoil-like arrangement)
            ang = [pi/2, 7*pi/6, 11*pi/6];
            cxy = P.rc * [cos(ang(:)), sin(ang(:))];

            % Internal equipotential regions (conductors)
            % NOTE: labels A/B/C can be re-mapped as needed for your figures.
            A = ((X - cxy(2,1)).^2 + (Y - cxy(2,2)).^2) <= P.rcond^2;
            B = ((X - cxy(1,1)).^2 + (Y - cxy(1,2)).^2) <= P.rcond^2;
            C = ((X - cxy(3,1)).^2 + (Y - cxy(3,2)).^2) <= P.rcond^2;

            chi(outer | A | B | C) = true;
            V(outer) = P.V0; % grounded screen
            V(A) = P.VA; V(B) = P.VB; V(C) = P.VC;

            plotNaNMask = outer; % hide outside of cable cross-section in plots

        case 2
            % Trench: grounded rectangular boundary (all edges)
            chi(1,:) = true; chi(end,:) = true; chi(:,1) = true; chi(:,end) = true;
            V(chi) = P.V0;

            % Conductors placed along a horizontal row at y = B - C
            yc = P.B - P.C;
            xc = [P.A/2 - P.D, P.A/2, P.A/2 + P.D];

            A = ((X - xc(1)).^2 + (Y - yc).^2) <= P.rcond^2;
            B = ((X - xc(2)).^2 + (Y - yc).^2) <= P.rcond^2;
            C = ((X - xc(3)).^2 + (Y - yc).^2) <= P.rcond^2;

            chi(A | B | C) = true;
            V(A) = P.VA; V(B) = P.VB; V(C) = P.VC;

        case 3
            % Rectangular box with Dirichlet boundary values
            chi(1,:) = true; chi(end,:) = true; chi(:,1) = true; chi(:,end) = true;

            if isfield(P,'boxConst') && P.boxConst
                V(chi) = P.V0;
            else
                % Assign each side independently (still Dirichlet)
                V(1,:)   = P.Vbottom;
                V(end,:) = P.Vtop;
                V(:,1)   = P.Vleft;
                V(:,end) = P.Vright;
            end

        case 4
            % Coaxial-like: grounded outer circle plus three internal conductors
            r = hypot(X, Y);
            outer = (r >= P.Rout);

            ang = [pi/2, 7*pi/6, 11*pi/6];
            cxy = P.rcent * [cos(ang(:)), sin(ang(:))];

            B = ((X - cxy(1,1)).^2 + (Y - cxy(1,2)).^2) <= P.rcond^2;
            A = ((X - cxy(2,1)).^2 + (Y - cxy(2,2)).^2) <= P.rcond^2;
            C = ((X - cxy(3,1)).^2 + (Y - cxy(3,2)).^2) <= P.rcond^2;

            chi(outer | A | B | C) = true;
            V(outer) = P.V0;
            V(A) = P.VA; V(B) = P.VB; V(C) = P.VC;

            plotNaNMask = outer;
    end
end

function [V, iters, dmax, hist] = laplaceFixedPoint(V, chi, hx, hy, epsV, kmax)
% Solves the discrete Laplace equation by fixed-point iterations on free nodes:
%   - Updates are computed from a 5-point finite-difference stencil
%     (paper: Eq. \eqref{eq:update_general}, derived from Eq. \eqref{eq:fdm_anisotropic})
%   - Fixed nodes (chi=1) keep their prescribed Dirichlet values
%     (paper: Eq. \eqref{eq:mask})
%   - Stop when ||ΔV||_inf < epsV over free nodes
%     (paper: Eq. \eqref{eq:stop_delta})
%
% Implementation detail:
%   This is a Jacobi-style update: V^{k+1} is built from V^{k} (stored as Vold).

    hx2 = hx*hx; hy2 = hy*hy;
    denom = 2*(hx2 + hy2); % denominator in Eq. \eqref{eq:update_general}

    freeAll = ~chi;
    hist = zeros(kmax,1);

    for k = 1:kmax
        Vold = V;
        Vnew = Vold;

        % Compute the finite-difference update V* on interior nodes (excluding edges)
        % Eq. \eqref{eq:update_general}:
        % V*_{i,j} = [ (V_{i+1,j}+V_{i-1,j}) hy^2 + (V_{i,j+1}+V_{i,j-1}) hx^2 ] / [ 2(hx^2+hy^2) ]
        Vstar = ((Vold(2:end-1,3:end) + Vold(2:end-1,1:end-2))*hy2 + ...
                 (Vold(3:end,2:end-1) + Vold(1:end-2,2:end-1))*hx2) / denom;

        % Apply updates only where nodes are free (chi=0)
        freeInterior = ~chi(2:end-1,2:end-1);
        tmp = Vnew(2:end-1,2:end-1);
        tmp(freeInterior) = Vstar(freeInterior);
        Vnew(2:end-1,2:end-1) = tmp;

        % Enforce Dirichlet values on fixed nodes exactly (chi=1)
        Vnew(chi) = Vold(chi);

        % Convergence check: ||ΔV||_inf on free nodes (paper: Eq. \eqref{eq:stop_delta})
        dV = abs(Vnew(freeAll) - Vold(freeAll));
        dmax = max(dV, [], 'all');
        hist(k) = dmax;

        V = Vnew;

        if dmax < epsV
            iters = k;
            hist = hist(1:k);
            return;
        end
    end

    iters = kmax;
end

function [Ex, Ey] = electricField(V, hx, hy)
% Electric field post-processing:
%   E = -∇V (paper: Eq. \eqref{eq:E_def})
% In a uniform grid, Ex and Ey correspond to finite differences of V.
% MATLAB gradient uses spacing (hy, hx) for (rows, cols), returning:
%   dVdy ~ ∂V/∂y and dVdx ~ ∂V/∂x
% Then:
%   Ex = -dVdx , Ey = -dVdy
% (paper: Eq. \eqref{eq:Ex_disc}--\eqref{eq:Ey_disc}, central differences in the interior)
    [dVdy, dVdx] = gradient(V, hy, hx);
    Ex = -dVdx;
    Ey = -dVdy;
end

function cmap = getCmap()
% Chooses a colormap for plotting. Uses turbo if available, else parula.
    if exist('turbo','file') == 2
        cmap = turbo;
    else
        cmap = parula;
    end
end