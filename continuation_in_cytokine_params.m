function heatmap_nu1
% HEATMAP_NU1
%
% Build a 2D map in (IFNg, IL-4) of the number of distinct attracting
% steady states in your 7D non-dimensional macrophage model
% for Hill exponent nu = 1 (non-ultrasensitive cybernetic control).
%
% x = [cat; inos; arg1; arg_int; NO; orn; arg_ext]

clear; clc;

%% 1. Global settings -----------------------------------------------------

params.epsd = 1e-8;
params.h    = 3.0;   % Hill exponent nu = 1 (no ultrasensitivity)

% Moderate / physiologically-inspired parameter set
% (these are the ones you used in your "Processes" draft code)
params.alpha_cat   = 0.01;
params.alpha_inos  = 0.01;
params.alpha_arg1  = 0.01;

params.beta_cat    = 0.2;
params.beta_inos   = 0.2;
params.beta_arg1   = 0.2;

params.kE_cat      = 1;
params.kE_inos     = 1;
params.kE_arg1     = 1;

params.KE_cat      = 10;
params.KE_inos     = 10;
params.KE_arg1     = 10;

params.mu_cat        = 1.5;
params.K_cat         = 10;
params.gamma_arg_int = 0.02;

params.mu_inos       = 1;
params.K_inos        = 10;
params.gamma_NO      = 0.05;

params.mu_arg1       = 1;
params.K_arg1        = 10;
params.gamma_orn     = 0.02;

params.K_perf        = 0.02;
params.Ab            = 100;
params.gamma_arg_ext = 0.02;

% Time scaling as in your non-dimensionalization
params.tau = 1/params.beta_cat;

% Derive normalization constants for non-dimensional state variables
params = compute_normalization_constants(params);

%% 2. Cytokine grids ------------------------------------------------------

% You can increase these if it runs fast enough
ifng_vals = linspace(0,1,15);   % IFNγ non-dim
il4_vals  = linspace(0,1,15);   % IL-4 non-dim

nI = numel(ifng_vals);
nL = numel(il4_vals);

% num_attr(i,j) = number of distinct attractors at (IFNg_i, IL4_j)
num_attr = zeros(nI,nL);

%% 3. Initial conditions for multi-attractor probing ----------------------

% A small, structured set of initial conditions: 
% near-zero, M1-like, M2-like, mid.
% x = [cat, inos, arg1, arg_int, NO, orn, arg_ext]
ICs = [
    1e-4  1e-4  1e-4  1e-4  1e-4  1e-4  1.0;  % almost zero
    1.0   0.9   0.1   0.8   0.8   0.2   1.0;  % M1-ish
    1.0   0.1   0.9   0.8   0.2   0.8   1.0;  % M2-ish
    0.5   0.5   0.5   0.5   0.5   0.5   1.0   % mid state
];
nIC = size(ICs,1);

%% 4. ODE solver settings -------------------------------------------------

T_end   = 300;  % final non-dimensional time (moderate, not insane)
odeOpts = odeset('RelTol',1e-5, 'AbsTol',1e-7);

tol_ss      = 1e-5;  % ||F(x_end)|| threshold to call it a steady state
tol_cluster = 1e-3;  % distance threshold to consider two states "same attractor"

%% 5. Main loop over cytokine grid ---------------------------------------

for i = 1:nI
    for j = 1:nL
        p_ifng = ifng_vals(i);
        p_il4  = il4_vals(j);

        attractors = [];  % columns = distinct attractors found here

        for s = 1:nIC
            x0 = ICs(s,:)';

            % Integrate non-dimensional ODE forward in time
            rhs = @(t,x) F_param_nd(x, p_ifng, p_il4, params);
            try
                [~, X] = ode15s(rhs, [0 T_end], x0, odeOpts);
            catch
                % solver failed (rare, but possible); skip this IC
                continue;
            end

            x_end = X(end,:)';
            x_end = max(x_end, 0);  % enforce non-negativity

            % Check if this looks like a steady state
            F_end = F_param_nd(x_end, p_ifng, p_il4, params);
            if norm(F_end, 2) > tol_ss
                % not close enough to a fixed point; ignore
                continue;
            end

            % Check if this steady state is new (cluster in state space)
            if isempty(attractors)
                attractors = x_end;
            else
                d = sqrt(sum( (attractors - x_end).^2, 1 ));
                if all(d > tol_cluster)
                    attractors(:,end+1) = x_end; %#ok<AGROW>
                end
            end
        end

        num_attr(i,j) = size(attractors,2);
    end
end

%% 6. Plot heatmap --------------------------------------------------------

figure;
imagesc(il4_vals, ifng_vals, num_attr);
set(gca,'YDir','normal');
xlabel('IL-4 (non-dim)');
ylabel('IFN\gamma (non-dim)');
title(sprintf('Number of attracting steady states (\\nu = %.1f)', params.h));

% Discrete colorbar 0–3
caxis([0 3]);
cb = colorbar;
cb.Ticks = 0:3;
cb.TickLabels = {'0','1','2','3'};
colormap(parula(4));

end % ---- end main function ----


%% =======================================================================
function params = compute_normalization_constants(params)
% Compute non-dimensional max scales for each state, in the same spirit
% as your previous code. These set the scaling of "x" in the model.

cat_0     = 0.001;
inos_0    = 0.001;
arg1_0    = 0.001;
arg_int_0 = 0.001;
NO_0      = 0.001;
orn_0     = 0.001;
arg_ext_0 = 1.0;

cat_max = max(cat_0, ...
    (params.alpha_cat + 2*(params.kE_cat/(params.KE_cat+1))) / params.beta_cat);

inos_max = max(inos_0, ...
    (params.alpha_inos + params.kE_inos/(params.KE_inos+1)) / params.beta_inos);

arg1_max = max(arg1_0, ...
    (params.alpha_arg1 + params.kE_arg1/(params.KE_arg1+1)) / params.beta_arg1);

arg_ext_max = max(arg_ext_0, ...
    (params.K_perf * params.Ab) / (params.K_perf + params.gamma_arg_ext));

arg_int_max = max(arg_int_0, ...
    (params.mu_cat * cat_max * arg_ext_max) / ...
    (params.gamma_arg_int * (params.K_cat + arg_ext_max)));

NO_max = max(NO_0, ...
    (params.mu_inos * inos_max * arg_int_max) / ...
    (params.gamma_NO * (params.K_inos + arg_int_max)));

orn_max = max(orn_0, ...
    (params.mu_arg1 * arg1_max * arg_int_max) / ...
    (params.gamma_orn * (params.K_arg1 + arg_int_max)));

params.cat_max     = cat_max;
params.inos_max    = inos_max;
params.arg1_max    = arg1_max;
params.arg_ext_max = arg_ext_max;
params.arg_int_max = arg_int_max;
params.NO_max      = NO_max;
params.orn_max     = orn_max;
end


%% =======================================================================
function F = F_param_nd(x, p_ifng, p_il4, params)
% Non-dimensional vector field with constant IFNg, IL-4.
% x = [cat; inos; arg1; arg_int; NO; orn; arg_ext]

cat     = x(1);
inos    = x(2);
arg1    = x(3);
arg_int = x(4);
NO      = x(5);
orn     = x(6);
arg_ext = x(7);

cat_max     = params.cat_max;
inos_max    = params.inos_max;
arg1_max    = params.arg1_max;
arg_ext_max = params.arg_ext_max;
arg_int_max = params.arg_int_max;
NO_max      = params.NO_max;
orn_max     = params.orn_max;
h           = params.h;
epsd        = params.epsd;

% ---- cybernetic fluxes and weights (Ramki-style) -----------------------
% physical fluxes (dimensional up to factor)
r_inos = (inos .* inos_max .* params.mu_inos .* arg_int .* arg_int_max) ./ ...
         (params.K_inos + arg_int * arg_int_max);

r_arg1 = (arg1 .* arg1_max .* params.mu_arg1 .* arg_int .* arg_int_max) ./ ...
         (params.K_arg1 + arg_int * arg_int_max);

% apply Hill exponent nu = h (here h=1 => no change, but code supports h>1)
r_inos_h = r_inos.^h;
r_arg1_h = r_arg1.^h;

den_u  = max(epsd, r_inos_h + r_arg1_h);
u_inos = r_inos_h ./ den_u;
u_arg1 = r_arg1_h ./ den_u;

den_v  = max(epsd, max(r_inos_h, r_arg1_h));
v_inos = r_inos_h ./ den_v;
v_arg1 = r_arg1_h ./ den_v;

F = zeros(7,1);

% 1) cat
F(1) = (params.tau*params.alpha_cat)/cat_max + ...
       (params.tau/cat_max)*(((params.kE_cat*p_ifng)/(params.KE_cat+p_ifng)) + ...
                             ((params.kE_cat*p_il4)/(params.KE_cat+p_il4))) ...
       - cat;

% 2) inos
F(2) = (params.tau*params.alpha_inos)/inos_max + ...
       (params.tau/inos_max)*((params.kE_inos * p_ifng * u_inos) / ...
                              (params.KE_inos + p_ifng)) ...
       - params.beta_inos * inos * params.tau;

% 3) arg1
F(3) = (params.tau*params.alpha_arg1)/arg1_max + ...
       (params.tau/arg1_max)*((params.kE_arg1 * p_il4 * u_arg1) / ...
                              (params.KE_arg1 + p_il4)) ...
       - params.beta_arg1 * arg1 * params.tau;

% 4) arg_int
F(4) = ...
    (params.tau*cat_max*cat*params.mu_cat*arg_ext*arg_ext_max) ./ ...
        (arg_int_max*params.K_cat + arg_int_max*arg_int_max*arg_ext) ...
  - (params.tau*v_inos .* inos .* arg_int * inos_max * params.mu_inos * arg_int_max) ./ ...
        (arg_int_max*params.K_inos + arg_int_max*arg_int_max*arg_int) ...
  - (params.tau*v_arg1 .* arg1 .* arg_int * arg1_max * params.mu_arg1 * arg_int_max) ./ ...
        (arg_int_max*params.K_arg1 + arg_int_max*arg_int_max*arg_int) ...
  - params.tau*params.gamma_arg_int * arg_int;

% 5) NO
F(5) = ...
    (params.tau*v_inos .* inos .* arg_int * inos_max * params.mu_inos * arg_int_max) ./ ...
        (NO_max*params.K_inos + NO_max*arg_int_max*arg_int) ...
  - params.tau*params.gamma_NO * NO;

% 6) orn
F(6) = ...
    (params.tau*v_arg1 .* arg1 .* arg_int * arg1_max * params.mu_arg1 * arg_int_max) ./ ...
        (orn_max*params.K_arg1 + orn_max*arg_int_max*arg_int) ...
  - params.tau*params.gamma_orn * orn;

% 7) arg_ext
F(7) = ...
    (params.tau/arg_ext_max)*params.K_perf*(params.Ab - arg_ext*arg_ext_max) ...
  - (params.tau*cat_max*params.mu_cat*arg_ext_max*cat.*arg_ext) ./ ...
        (arg_ext_max*params.K_cat + arg_ext_max*arg_ext_max*arg_ext) ...
  - params.gamma_arg_ext * params.tau * arg_ext;

end
