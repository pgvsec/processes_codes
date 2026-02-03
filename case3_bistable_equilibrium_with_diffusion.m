%% Population-balance simulation with KDE-based marginals (CORRECTED + drift_fix)
clear; clc;

% x vector = [cat inos arg1 arg_int NO orn]

% time-step size (dimensional)
h  = 0.01;
tf = 196; % final time in hours

% general parameters
params.epsd        = 1e-8;
params.slow_factor = 0.75;

% intializing vectors for a single-cell ODE (used for scales)
t        = zeros(tf/h,1);       t(1)       = 0; % time
cat      = zeros(tf/h,1);       cat(1)     = 0.001;
inos     = zeros(tf/h,1);       inos(1)    = 0.001;
arg1     = zeros(tf/h,1);       arg1(1)    = 0.001;
arg_int  = zeros(tf/h,1);       arg_int(1) = 0.001;
NO       = zeros(tf/h,1);       NO(1)      = 0.001;
orn      = zeros(tf/h,1);       orn(1)     = 0.001;
arg_ext  = zeros(tf/h,1);       arg_ext(1) = 1;

% cytokine signals (just constants here; could be gaussians in time)
ifng_peak   = 24;   %#ok<NASGU> % in hours
il4_peak    = 24*3; %#ok<NASGU>
ifng_spread = 10;   %#ok<NASGU>
il4_spread  = 30;   %#ok<NASGU>
ifng        = .6;   %#ok<NASGU>
il4         = .0;   %#ok<NASGU>

%% Pack enzyme parameters into struct

% basal synthesis rates (nM·h⁻¹)
params.alpha.cat   = 0.01;
params.alpha.inos  = 0.01;
params.alpha.arg1  = 0.01;

% first-order degradation rates (h⁻¹)
params.beta.cat    = 0.2;
params.beta.inos   = 0.2;
params.beta.arg1   = 0.2;

% maximal inducible synthesis (nM·h⁻¹)
params.kE.cat      = 1;
params.kE.inos     = 1;
params.kE.arg1     = 1;

% EC50-like parameters (nM)
params.KE.cat      = 10;
params.KE.inos     = 10;
params.KE.arg1     = 10;

%% Metabolite/cytokine parameters

params.mu.cat        = 1.5;
params.K.cat         = 10;
params.gamma.arg_int = 0.2;

params.mu.inos       = 1;
params.K.inos        = 10;
params.gamma.NO      = 0.5;

params.mu.arg1       = 1;
params.K.arg1        = 10;
params.gamma.orn     = 0.5;

params.K_perf        = 0.02;
params.Ab            = 100; % blood levels in micro-molar
params.gamma.arg_ext = 0.02;

params.tau           = 1/params.beta.cat; % time scale

%% Non-dimensional cytokine levels (here taken as constants)
ifng_nondim = 0.8;
il4_nondim  = 0.8;

%% Compute max scales for non-dimensionalization (rough upper bounds)

cat_max     = max(cat(1), (params.alpha.cat+2*(params.kE.cat/(params.KE.cat+1)))/params.beta.cat);
inos_max    = max(inos(1), (params.alpha.inos+((params.kE.inos)/(params.KE.inos+1)))/params.beta.inos);
arg1_max    = max(arg1(1), (params.alpha.arg1+((params.kE.arg1)/(params.KE.arg1+1)))/params.beta.arg1);
arg_ext_max = max(arg_ext(1), (params.K_perf*params.Ab)/(params.K_perf+params.gamma.arg_ext));
arg_int_max = max(arg_int(1), (params.mu.cat*cat_max*arg_ext_max)/(params.gamma.arg_int*(params.K.cat+arg_ext_max)));
NO_max      = max(NO(1), (params.mu.inos*inos_max*arg_int_max)/(params.gamma.NO*(params.K.inos+arg_int_max)));
orn_max     = max(orn(1),(params.mu.arg1*arg1_max*arg_int_max)/(params.gamma.orn*(params.K.arg1+arg_int_max)));

%% Cybernetic weights (non-dimensional, vectorized form)

nu = 2;  % Hill exponent

u_inos = @(cat_nd,inos_nd,arg1_nd,arg_int_nd,NO_nd,orn_nd) ...
    (((inos_nd.*inos_max.*params.mu.inos.*arg_int_nd*arg_int_max)./(params.K.inos+arg_int_nd*arg_int_max)).^nu) ...
    ./ max(params.epsd, ...
        ((inos_nd.*inos_max.*params.mu.inos.*arg_int_nd*arg_int_max)./(params.K.inos+arg_int_nd*arg_int_max)).^nu + ...
        ((arg1_nd.*arg1_max.*params.mu.arg1.*arg_int_nd*arg_int_max)./(params.K.arg1+arg_int_nd*arg_int_max)).^nu);

u_arg1 = @(cat_nd,inos_nd,arg1_nd,arg_int_nd,NO_nd,orn_nd) ...
    (((arg1_nd.*arg1_max.*params.mu.arg1.*arg_int_nd*arg_int_max)./(params.K.arg1+arg_int_nd*arg_int_max)).^nu) ...
    ./ max(params.epsd, ...
        ((inos_nd.*inos_max.*params.mu.inos.*arg_int_nd*arg_int_max)./(params.K.inos+arg_int_nd*arg_int_max)).^nu + ...
        ((arg1_nd.*arg1_max.*params.mu.arg1.*arg_int_nd*arg_int_max)./(params.K.arg1+arg_int_nd*arg_int_max)).^nu);

v_inos = @(cat_nd,inos_nd,arg1_nd,arg_int_nd,NO_nd,orn_nd) ...
    (((inos_nd.*inos_max.*params.mu.inos.*arg_int_nd*arg_int_max)./(params.K.inos+arg_int_nd*arg_int_max)).^nu) ...
    ./ max(params.epsd, max( ...
        ((inos_nd.*inos_max.*params.mu.inos.*arg_int_nd*arg_int_max)./(params.K.inos+arg_int_nd*arg_int_max)).^nu, ...
        ((arg1_nd.*arg1_max.*params.mu.arg1.*arg_int_nd*arg_int_max)./(params.K.arg1+arg_int_nd*arg_int_max)).^nu));

v_arg1 = @(cat_nd,inos_nd,arg1_nd,arg_int_nd,NO_nd,orn_nd) ...
    ((arg1_nd.*arg1_max.*params.mu.arg1.*arg_int_nd*arg_int_max)./(params.K.arg1+arg_int_nd*arg_int_max)).^nu ...
    ./ max(params.epsd, max( ...
        ((inos_nd.*inos_max.*params.mu.inos.*arg_int_nd*arg_int_max)./(params.K.inos+arg_int_nd*arg_int_max)).^nu, ...
        ((arg1_nd.*arg1_max.*params.mu.arg1.*arg_int_nd*arg_int_max)./(params.K.arg1+arg_int_nd*arg_int_max)).^nu));

% Hill repressions (currently turned off = 1)
Ho_fun = @(P) 1; % Orn → iNOS repression
Hn_fun = @(P) 1; % NO  → Arg1 repression

%% --- Vectorized drift for particles (non-dimensional) ---
% F_particles_vec: returns an array that might be N x 6 or flattened.
% We will *always* pass its output through drift_fix(...) before using.
F_particles_vec = @(t_nd, P, z_scalar) params.slow_factor * [ ...
    % F1: cat
    (params.tau*params.alpha.cat)/cat_max ...
      + (params.tau/cat_max)*(((params.kE.cat*ifng_nondim)/(params.KE.cat+ifng_nondim)) ...
                             +((params.kE.cat*il4_nondim )./(params.KE.cat+il4_nondim ))) ...
      - P(:,1), ...
    % F2: iNOS
    (params.tau*params.alpha.inos)/inos_max ...
      + Ho_fun(P).*(params.tau/inos_max).* ...
        ((params.kE.inos*ifng_nondim .* u_inos(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)))./(params.KE.inos+ifng_nondim)) ...
      - params.beta.inos*P(:,2)*params.tau, ...
    % F3: Arg1
    (params.tau*params.alpha.arg1)/arg1_max ...
      + Hn_fun(P).*(params.tau/arg1_max).* ...
        ((params.kE.arg1*il4_nondim .* u_arg1(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)))./(params.KE.arg1+il4_nondim)) ...
      - params.beta.arg1*P(:,3)*params.tau, ...
    % F4: Arg_int
    (params.tau*cat_max*P(:,1)*params.mu.cat*z_scalar*arg_ext_max)./(arg_int_max*params.K.cat + arg_int_max*arg_int_max*z_scalar) ...
      - (params.tau*v_inos(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,2).*P(:,4)*inos_max*params.mu.inos*arg_int_max) ...
        ./(arg_int_max*params.K.inos + arg_int_max*arg_int_max*P(:,4)) ...
      - (params.tau*v_arg1(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,3).*P(:,4)*arg1_max*params.mu.arg1*arg_int_max) ...
        ./(arg_int_max*params.K.arg1 + arg_int_max*arg_int_max*P(:,4)) ...
      - params.tau*params.gamma.arg_int*P(:,4), ...
    % F5: NO
    (params.tau*v_inos(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,2).*P(:,4)*inos_max*params.mu.inos*arg_int_max) ...
        ./(NO_max*params.K.inos + NO_max*arg_int_max*P(:,4)) ...
      - params.tau*params.gamma.NO*P(:,5), ...
    % F6: Orn
    (params.tau*v_arg1(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,3).*P(:,4)*arg1_max*params.mu.arg1*arg_int_max) ...
        ./(orn_max*params.K.arg1 + orn_max*arg_int_max*P(:,4)) ...
      - params.tau*params.gamma.orn*P(:,6) ...
];

%% Initial particle cloud X0 ~ truncated multivariate gaussian in [0,1]^6

% x vector = [cat inos arg1 arg_int NO orn]
mu0 = 50*[0.01 0.01 0.01 0.01 0.01 0.01]; % initial mean (can adjust)

alpha = 0.5;
sigma = zeros(1,6);
for i = 1:length(mu0)
    sigma(i) = alpha*min(mu0(i),1-mu0(i));
end

% correlation matrix R
R = eye(6); shrink = 0.7;

% correlations
R(1,4) = -0.5*shrink; R(4,1) = -0.5*shrink;  % CAT <-> Arg_int
R(2,3) = -0.5*shrink; R(3,2) = -0.5*shrink;  % iNOS <-> Arg1
R(2,4) = -0.5*shrink; R(4,2) = -0.5*shrink;  % iNOS <-> Arg_int
R(3,4) = -0.5*shrink; R(4,3) = -0.5*shrink;  % Arg1 <-> Arg_int
R(2,5) =  0.5*shrink; R(5,2) =  0.5*shrink;  % iNOS <-> NO
R(3,6) =  0.5*shrink; R(6,3) =  0.5*shrink;  % Arg1 <-> Orn

D     = diag(sigma);
Sigma = D * R * D;

% PS definiteness check
eigsSigma = eig((Sigma+Sigma')/2);
if min(eigsSigma) <= 0
    error('Covariance Sigma is not positive semidefinite. Adjust R/sigma.');
end

Np  = 1e4;     % number of particles
dim = 6;
X0  = zeros(Np, dim);
have = 0;

rng(1);
batch = max(1000, ceil(0.02*Np));

while have < Np
    draw  = mvnrnd(mu0, Sigma, batch);
    inBox = all(draw >= 0 & draw <= 1, 2);
    keep  = draw(inBox,:);
    nk    = min(size(keep,1), Np - have);
    if nk > 0
        X0(have+1:have+nk, :) = keep(1:nk, :);
        have = have + nk;
    end
end

P_current = max(0, min(1, X0)); % initial internal states in [0,1]^6
z_current = arg_ext(1);         % initial external arginine (non-dim here)

%% Initial weights & volumes

% Choose total mass M0 arbitrarily (e.g. 1)
M0    = 1;
w_vec = (M0 / Np) * ones(Np, 1);   % each particle gets equal initial mass
V_vec = (1 / Np) * ones(Np, 1);    % volumes normalized
V_current = V_vec;

%% Source term inputs

x_b       = [0 0 0 0 0 0];        % baseline state for source
sigma_hat = 0.04*ones(1,6);       % Gaussian widths
A_fun     = @(t_dim) (0.02 + 0.01*sin(0.1*t_dim)); % dimensional amplitude

%% APM options

APM_opts.N_target     = Np;
APM_opts.merge_factor =  0;%0.5;
APM_opts.split_factor = 2.0;
APM_opts.d_max        = 0.02;
APM_opts.split_jitter = 0.05;

%% aux struct for analytic divergence

aux.inos_max    = inos_max;
aux.arg1_max    = arg1_max;
aux.arg_int_max = arg_int_max;
aux.ifng_nondim = ifng_nondim;
aux.il4_nondim  = il4_nondim;
aux.kappa       = 20;
aux.nu          = nu;

%% Diffusion kernel parameters

params.D_int   = 0.01*params.slow_factor; % internal diffusion strength
params.eps_kde = 0.15;                    % kernel bandwidth in [0,1]
params.u_floor = 1e-3;                    % floor on n_eps
params.K_NO = 0.4;  % initial guess; will be updated during simulation

%% External arginine evolution (scalar)

R_ext_handle       = @(z_scalar) ((params.K_perf*params.tau)/arg_ext_max)*(params.Ab - z_scalar*arg_ext_max) - params.gamma.arg_ext*z_scalar*params.tau;
uptake_rate_handle = @(P, z_scalar) (params.tau*P(:,1)*cat_max*params.mu.cat*arg_ext_max*z_scalar) ...
                                  ./ (arg_int_max*params.K.cat+arg_int_max*arg_ext_max*z_scalar);

F_external_z = @(t_nd, P, z_scalar, w_vec) ...
    -sum(w_vec(:) .* uptake_rate_handle(P, z_scalar)) + R_ext_handle(z_scalar);

%% Time integration (non-dimensional time)

tn        = 0;
dt        = 0.05;
tn_final  = tf / params.tau;
N_steps   = ceil(tn_final / dt);

t_save = [1, 4, 20, 40, 80, 120, 160, 190] / params.tau;
Nt     = numel(t_save);

P_snap = cell(1,Nt);
w_snap = cell(1,Nt);
V_snap = cell(1,Nt);
z_snap = cell(1,Nt);
t_snap = nan(1,Nt);
saved  = false(1,Nt);

inos_x  = cell(1,Nt);
inos_f  = cell(1,Nt);
arg1_x  = cell(1,Nt);
arg1_f  = cell(1,Nt);

cloud_inos_arg1 = cell(1,Nt);
cloud_w         = cell(1,Nt);

inos_mean = nan(1,Nt);
arg1_mean = nan(1,Nt);

% save t=0 snapshot if requested
idx0 = find(abs(t_save - 0) < 1e-12, 1);
if ~isempty(idx0)
    P_snap{idx0} = P_current;
    w_snap{idx0} = w_vec(:);
    V_snap{idx0} = V_current(:);
    z_snap{idx0} = z_current;
    t_snap(idx0) = 0;
    saved(idx0)  = true;
end

%% Main RK3 loop
for n = 1:N_steps

    % sanity: keep particles and z within valid ranges
    P_current = min(1, max(0, P_current));
    z_current = max(0, z_current);

    % (optional sanity check)
    if size(P_current,1) ~= numel(w_vec)
        error('Mismatch: size(P_current,1)=%d, numel(w_vec)=%d', ...
              size(P_current,1), numel(w_vec));
    end

    % ==================== STAGE 1 ====================
    F_drift_n_raw = F_particles_vec(tn, P_current, z_current);
    F_drift_n     = drift_fix(F_drift_n_raw, P_current);          % force N x 6
    A_diff_n      = diffusion_drift_internal_vec(P_current, w_vec, params);
    kP_n          = F_drift_n + A_diff_n;

    kV_n = analytic_divergence(tn, P_current, z_current, params, aux).* V_current;
    kw_n = -(DeathTerm(P_current, params, NO_max).*w_vec) ...
           + SourceTerm(P_current, x_b, sigma_hat, A_fun, tn*params.tau, params).*V_current;
    kz_n = F_external_z(tn, P_current, z_current, w_vec);

    P1 = P_current + dt * kP_n;
    V1 = V_current + dt * kV_n;
    w1 = w_vec     + dt * kw_n;  w1 = max(0, w1(:));
    z1 = z_current + dt * kz_n;

    % ==================== STAGE 2 ====================
    tn_plus_dt = tn + dt;

    F_drift_1_raw = F_particles_vec(tn_plus_dt, P1, z1);
    F_drift_1     = drift_fix(F_drift_1_raw, P1);
    A_diff_1      = diffusion_drift_internal_vec(P1, w1, params);
    kP_1          = F_drift_1 + A_diff_1;

    kV_1 = analytic_divergence(tn_plus_dt, P1, z1, params, aux).* V1;
    kw_1 = -(DeathTerm(P1, params, NO_max).*w1) ...
           + SourceTerm(P1, x_b, sigma_hat, A_fun, tn_plus_dt*params.tau, params).*V1;
    kz_1 = F_external_z(tn_plus_dt, P1, z1, w1);

    P2 = (3/4)*P_current + (1/4)*P1 + (1/4)*dt * kP_1;
    V2 = (3/4)*V_current + (1/4)*V1 + (1/4)*dt * kV_1;
    w2 = (3/4)*w_vec     + (1/4)*w1 + (1/4)*dt * kw_1; w2 = max(0, w2(:));
    z2 = (3/4)*z_current + (1/4)*z1 + (1/4)*dt * kz_1;

    % ==================== STAGE 3 ====================
    F_drift_2_raw = F_particles_vec(tn_plus_dt, P2, z2);
    F_drift_2     = drift_fix(F_drift_2_raw, P2);
    A_diff_2      = diffusion_drift_internal_vec(P2, w2, params);
    kP_2          = F_drift_2 + A_diff_2;

    kV_2 = analytic_divergence(tn_plus_dt, P2, z2, params, aux).* V2;
    kw_2 = -(DeathTerm(P2, params, NO_max).*w2) ...
           + SourceTerm(P2, x_b, sigma_hat, A_fun, tn_plus_dt*params.tau, params).*V2;
    kz_2 = F_external_z(tn_plus_dt, P2, z2, w2);

    P_next = (1/3)*P_current + (2/3)*P2 + (2/3)*dt * kP_2;
    V_next = (1/3)*V_current + (2/3)*V2 + (2/3)*dt * kV_2;
    w_next = (1/3)*w_vec     + (2/3)*w2 + (2/3)*dt * kw_2;
    z_next = (1/3)*z_current + (2/3)*z2 + (2/3)*dt * kz_2;

    % update and clip
    P_current = min(1, max(0, P_next));
    V_current = max(0, V_next);
    w_vec     = max(0, w_next(:));
    z_current = max(0, z_next);

    tn = tn_plus_dt;

    % ========= SNAPSHOTS (particles, weights, z) =========
    tol_t  = 0.5*dt;
    idx_hit = find(~saved & abs(tn - t_save) <= tol_t, 1, 'first');
    if ~isempty(idx_hit)
        k_idx = idx_hit;
        P_snap{k_idx} = P_current;
        w_snap{k_idx} = w_vec(:);
        V_snap{k_idx} = V_current(:);
        z_snap{k_idx} = z_current;
        t_snap(k_idx) = tn;
        saved(k_idx)  = true;
    end

    % ========= KDE snapshots for iNOS and Arg1 =========
   % ------------------KDE snapshots 1D --------------------
tol_t  = 0.5*dt; % tolerance so tn 'matches' a requested save time
idx_hit = find(abs(tn - t_save) <= tol_t, 1, 'first');

if ~isempty(idx_hit)
    k = idx_hit;

    % total mass (measure)
    M = sum(w_vec);
    if M <= 0
        fprintf('t = %.2f (nd): total mass <= 0, skipping KDE.\n', tn);
        return;  % or: continue; inside your loop
    end

    %=============
        % ===== Basin-mass + NO diagnostics (REAL asymmetry) =====
    theta = 0.55;  % choose a midline once (0.5–0.6 typically)
    inos_now = P_current(:,2);
    NO_dim   = NO_max .* P_current(:,5);

    M_low  = sum(w_vec(inos_now < theta));
    M_high = sum(w_vec(inos_now > theta));

    % avoid empty sets
    if any(inos_now < theta)
        NO_low_mean = mean(NO_dim(inos_now < theta));
    else
        NO_low_mean = NaN;
    end
    if any(inos_now > theta)
        NO_high_mean = mean(NO_dim(inos_now > theta));
    else
        NO_high_mean = NaN;
    end
    % Update K_NO to sit between basin NO levels (only if both exist)
    if isfinite(NO_low_mean) && isfinite(NO_high_mean)
        params.K_NO = 0.5*(NO_low_mean + NO_high_mean);
    end

    fprintf('t=%.2f nd: M_low=%.3e, M_high=%.3e, ratio(high/low)=%.3f | NO_low=%.3g, NO_high=%.3g\n', ...
        tn, M_low, M_high, M_high/max(M_low,1e-30), NO_low_mean, NO_high_mean);
%================

    % ---------- pull coordinates ----------
    inos_vals = real(P_current(:,2));  % iNOS
    arg1_vals = real(P_current(:,3));  % Arg1

    % ---------- build base normalized weights ----------
    w_pdf = w_vec(:);        % mass weights
    w_pdf = w_pdf / M;       % now sum(w_pdf) = 1 (if no NaNs)

    % ---------- remove NaN/Inf / out-of-support points ----------
    % We require data in [0,1] and finite, and weights finite.
    valid_inos = isfinite(inos_vals) & isfinite(w_pdf);
    valid_arg1 = isfinite(arg1_vals) & isfinite(w_pdf);

    % (you can choose to require both coords valid, but here we do 1D cleanly)
    inos_vals = inos_vals(valid_inos);
    w_inos    = w_pdf(valid_inos);

    arg1_vals = arg1_vals(valid_arg1);
    w_arg1    = w_pdf(valid_arg1);

    % Clip strictly to [0,1] AFTER filtering
    inos_vals = min(1, max(0, inos_vals));
    arg1_vals = min(1, max(0, arg1_vals));

    % Remove any possible tiny roundoff violations again
    valid_inos2 = (inos_vals >= 0) & (inos_vals <= 1) & isfinite(inos_vals) & isfinite(w_inos);
    valid_arg12 = (arg1_vals >= 0) & (arg1_vals <= 1) & isfinite(arg1_vals) & isfinite(w_arg1);

    inos_vals = inos_vals(valid_inos2);
    w_inos    = w_inos(valid_inos2);
    arg1_vals = arg1_vals(valid_arg12);
    w_arg1    = w_arg1(valid_arg12);

    % If nothing left, skip KDE for this time
    if numel(inos_vals) < 5 || sum(w_inos) <= 0
        fprintf('t = %.2f (nd): not enough valid iNOS data for KDE, skipping.\n', tn);
    else
        % Renormalize weights for iNOS KDE
        w_inos = w_inos / sum(w_inos);

        fprintf('t = %.2f (nd), iNOS: Nvalid = %d, min = %.3e, max = %.3e\n', ...
                tn, numel(inos_vals), min(inos_vals), max(inos_vals));
        fprintf('   sum(w_iNOS) = %.6e\n', sum(w_inos));

        % evaluation grid restricted to [0,1]
xi_inos = linspace(0,1,200);

[f_inos_pdf, x_inos] = ksdensity(inos_vals, xi_inos, ...
    'Weights', w_inos, ...
    'BoundaryCorrection', 'reflection');

        fprintf('   integral pdf_iNOS ~ %.3e\n', trapz(x_inos, f_inos_pdf));

        % convert to number density: integrate ≈ total mass M
        inos_x{k} = x_inos;
        inos_f{k} = f_inos_pdf * M;

        dx_inos          = x_inos(2) - x_inos(1);
        inos_mean(k)     = sum(x_inos(:) .* f_inos_pdf(:)) * dx_inos;
    end

    if numel(arg1_vals) < 5 || sum(w_arg1) <= 0
        fprintf('t = %.2f (nd): not enough valid Arg1 data for KDE, skipping.\n', tn);
    else
        % Renormalize weights for Arg1 KDE
        w_arg1 = w_arg1 / sum(w_arg1);

        fprintf('t = %.2f (nd), Arg1: Nvalid = %d, min = %.3e, max = %.3e\n', ...
                tn, numel(arg1_vals), min(arg1_vals), max(arg1_vals));
        fprintf('   sum(w_Arg1) = %.6e\n', sum(w_arg1));

       xi_arg1 = linspace(0,1,200);

[f_arg1_pdf, x_arg1] = ksdensity(arg1_vals, xi_arg1, ...
    'Weights', w_arg1, ...
    'BoundaryCorrection', 'reflection');

        arg1_x{k} = x_arg1;
        arg1_f{k} = f_arg1_pdf * M;

        dx_arg1          = x_arg1(2) - x_arg1(1);
        arg1_mean(k)     = sum(x_arg1(:) .* f_arg1_pdf(:)) * dx_arg1;
    end

    % ----- subsampled cloud for visualization (optional) -----
    Nsamp    = min(2000, size(P_current,1));
    idx_samp = randperm(size(P_current,1), Nsamp);
    cloud_inos_arg1{k} = P_current(idx_samp, [2 3]);
    cloud_w{k}         = w_pdf(idx_samp);  % original normalized mass weights
end


    % ========= APM every 20 steps =========
     if mod(n,20) == 0
         [P_current, w_vec, V_current] = adaptive_particle_management(P_current, w_vec, V_current, APM_opts);
         P_current = min(1, max(0, P_current));
         w_vec     = max(0, w_vec(:));
         V_current = max(0, V_current(:));
     end

end % main loop

%% -------------------------------- Results and Visualizations ------------------------------

% 1) iNOS number-density snapshots
figure;
for k = 1:Nt
    if isempty(inos_x{k}), continue; end
    plot(inos_x{k}, inos_f{k}, 'DisplayName', sprintf('t = %.1f (nd)', t_save(k)));
    hold on;
end
xlabel('iNOS (nondimensional)');
ylabel('Number density (mass per unit state)');
title('iNOS number-density snapshots');
legend show; grid on;

% 1)* iNOS snapshots with inset figure
figure;

% -------- MAIN AXES: all times ----------
ax_main = axes;
hold(ax_main,'on');

for k = 1:Nt
    if isempty(inos_x{k}), continue; end
    plot(ax_main, inos_x{k}, inos_f{k}, ...
        'LineWidth',1.1,'DisplayName', sprintf('t = %.1f (nd)', t_save(k)));
end

xlabel(ax_main,'iNOS (nondimensional)');
ylabel(ax_main,'Number density (mass per unit state)');
title(ax_main,'iNOS number-density snapshots');
legend(ax_main,'show'); grid(ax_main,'on');
xlim(ax_main,[0 1]);

% -------- INSET: last 3 times only ----------
idx_last = max(1, Nt-2):Nt;   % indices for last 3 snapshots

% compute max value among last 3 to set vertical zoom
max_zoom = 0;
for k = idx_last
    if ~isempty(inos_f{k})
        max_zoom = max(max_zoom, max(inos_f{k}));
    end
end

ax_inset = axes('Position',[0.55 0.55 0.35 0.35]);  % [left bottom width height]
hold(ax_inset,'on');

for k = idx_last
    if isempty(inos_x{k}), continue; end
    plot(ax_inset, inos_x{k}, inos_f{k}, ...
         'LineWidth',1.1,'DisplayName', sprintf('t = %.1f (nd)', t_save(k)));
end

xlim(ax_inset,[0 1]);          % or e.g. [0.7 1] if you want to focus on M2 peak
ylim(ax_inset,[0 1.1*max_zoom]);  % vertical zoom so late peaks are tall
grid(ax_inset,'on');
box(ax_inset,'on');
set(ax_inset,'FontSize',8);    % smaller labels in inset
legend(ax_inset,'show','Location','northwest');


% 2) Arg1 number-density snapshots
figure;
for k = 1:Nt
    if isempty(arg1_x{k}), continue; end
    plot(arg1_x{k}, arg1_f{k}, 'DisplayName', sprintf('t = %.1f (nd)', t_save(k)));
    hold on;
end
xlabel('Arg1 (nondimensional)');
ylabel('Number density (mass per unit state)');
title('Arg1 number-density snapshots');
legend show; grid on;

% % 3) Particle cloud in (iNOS, Arg1) at different times
% figure;
% for k = 1:Nt
%     if isempty(cloud_inos_arg1{k}), continue; end
%     subplot(2, ceil(Nt/2), k);
%     pts = cloud_inos_arg1{k};
%     scatter(pts(:,1), pts(:,2), 8, cloud_w{k}, 'filled');
%     xlabel('iNOS (nd)'); ylabel('Arg1 (nd)');
%     title(sprintf('t = %.1f (nd)', t_save(k)));
%     xlim([0 1]); ylim([0 1]);
%     colorbar;
% end
% sgtitle('Particle cloud in (iNOS, Arg1) space');
idx_plot = [1 3 5 8];    % <-- choose any 4 indices from 1:Nt
n_plot   = numel(idx_plot);

figure;
for j = 1:n_plot
    k = idx_plot(j);                 % actual snapshot index

    if isempty(cloud_inos_arg1{k}), continue; end

    subplot(2,2,j);                  % always 2x2 layout
    pts = cloud_inos_arg1{k};
    scatter(pts(:,1), pts(:,2), 8, cloud_w{k}, 'filled');
    xlabel('iNOS (nd)');
    ylabel('Arg1 (nd)');
    title(sprintf('t = %.1f (nd)', t_save(k)));
    xlim([0 1]); ylim([0 1]);
    colorbar;
end
sgtitle('Particle cloud in (iNOS, Arg1) space');

% 4) Phase portrait at final time
t_phase = tn;
z_phase = z_current;

cat0 = mean(P_current(:,1));
q0   = mean(P_current(:,4));
NO0  = mean(P_current(:,5));
orn0 = mean(P_current(:,6));

[inos_grid, arg1_grid] = meshgrid(linspace(0,1,20), linspace(0,1,20));
Ng = numel(inos_grid);

P_grid = zeros(Ng, 6);
P_grid(:,1) = cat0;
P_grid(:,2) = inos_grid(:);
P_grid(:,3) = arg1_grid(:);
P_grid(:,4) = q0;
P_grid(:,5) = NO0;
P_grid(:,6) = orn0;

F_grid_raw = F_particles_vec(t_phase, P_grid, z_phase);
F_grid     = drift_fix(F_grid_raw, P_grid);

u_vec = reshape(F_grid(:,2), size(inos_grid));
v_vec = reshape(F_grid(:,3), size(arg1_grid));

figure;
quiver(inos_grid, arg1_grid, u_vec, v_vec); hold on;
pts_final = cloud_inos_arg1{end};
if ~isempty(pts_final)
    scatter(pts_final(:,1), pts_final(:,2), 10, 'r', 'filled');
end
xlabel('iNOS (nd)');
ylabel('Arg1 (nd)');
title(sprintf('Phase portrait in (iNOS, Arg1) at t = %.1f (nd)', t_phase));
xlim([0 1]); ylim([0 1]);
legend('vector field','particles');

% 5) Attractor candidate & Jacobian
w_norm_final = w_vec(:) / sum(w_vec);
x_star = sum(P_current .* w_norm_final, 1);
z_star = z_current;
t_star = t_phase;

fprintf('\n===== ATTRACTOR CANDIDATE x_star =====\n');
disp(x_star);

Fx_raw  = F_particles_vec(t_star, x_star, z_star);
Fx      = drift_fix(Fx_raw, x_star);
Fx_norm = norm(Fx);
fprintf('||F(x_star)|| = %.3e\n', Fx_norm);

J = jacobian_complex(@(t,X,z) drift_fix(F_particles_vec(t,X,z), X), x_star, t_star, z_star);
fprintf('\n===== EIGENVALUES OF JACOBIAN @ ATTRACTOR =====\n');
eigvals = eig(J);
disp(eigvals);

if all(real(eigvals) < 0)
    disp('=> The attractor is a LOCALLY STABLE sink.');
elseif any(real(eigvals) > 0)
    disp('=> The attractor is UNSTABLE or a saddle.');
else
    disp('=> The attractor has neutrally stable directions.');
end

%% --------------------------- Local Functions --------------------------------

function F = drift_fix(F_raw, P)
    % Make sure F has the same shape as P.
    if isvector(F_raw) && numel(F_raw) == numel(P)
        F = reshape(F_raw, size(P));   % e.g. 6000x1 -> 1000x6
    else
        F = F_raw;
    end
end

function div = analytic_divergence(t, P, z_scalar, params, aux) %#ok<INUSL>
    [N, dim] = size(P);
    if dim ~= 6
        error('analytic_divergence: P must be N x 6.');
    end

    c = P(:,1); %#ok<NASGU>
    i = P(:,2);
    a = P(:,3);
    q = P(:,4);
    n = P(:,5); %#ok<NASGU>
    o = P(:,6); %#ok<NASGU>

    tau    = params.tau;
    inos_max    = aux.inos_max;
    arg1_max    = aux.arg1_max;
    q_max       = aux.arg_int_max;
    mu_inos     = params.mu.inos;
    mu_arg1     = params.mu.arg1;
    K_inos      = params.K.inos;
    K_arg1      = params.K.arg1;
    beta_inos   = params.beta.inos;
    beta_arg1   = params.beta.arg1;
    gamma_q     = params.gamma.arg_int;
    gamma_NO    = params.gamma.NO;
    gamma_orn   = params.gamma.orn;
    kE_inos     = params.kE.inos;
    kE_arg1     = params.kE.arg1;
    KE_inos     = params.KE.inos;
    KE_arg1     = params.KE.arg1;
    kappa       = aux.kappa;
    nu          = aux.nu;
    ifng_val    = aux.ifng_nondim;
    il4_val     = aux.il4_nondim;
    tiny        = 1e-12;

    % dF1/dc
    dF1_dc = -ones(N,1);

    % Raw rates R1, R2
    A1   = i .* inos_max * mu_inos * q_max;
    den1 = K_inos + q * q_max;
    R1   = (A1 .* q) ./ (den1 + tiny);

    A2   = a .* arg1_max * mu_arg1 * q_max;
    den2 = K_arg1 + q * q_max;
    R2   = (A2 .* q) ./ (den2 + tiny);

    R1_q = A1 .* K_inos ./ (den1 + tiny).^2;
    R2_q = A2 .* K_arg1 ./ (den2 + tiny).^2;

    C1 = (inos_max * mu_inos * q .* q_max) ./ (den1 + tiny);
    C2 = (arg1_max * mu_arg1 * q .* q_max) ./ (den2 + tiny);

    % Hill-transformed rates
    R1h = R1.^nu;
    R2h = R2.^nu;
    S   = R1h + R2h + tiny;
    u_inos_loc = R1h ./ S; %#ok<NASGU>
    u_arg1_loc = R2h ./ S; %#ok<NASGU>

    dR1h_di = nu * R1.^(nu-1) .* C1;
    dR2h_da = nu * R2.^(nu-1) .* C2;

    du_inos_di = dR1h_di .* (R2h ./ (S.^2));
    du_arg1_da = dR2h_da .* (R1h ./ (S.^2));

    % iNOS and Arg1 diagonal derivatives
    B_inos = tau / inos_max * (kE_inos * ifng_val / (KE_inos + ifng_val + tiny));
    B_arg1 = tau / arg1_max * (kE_arg1 * il4_val / (KE_arg1 + il4_val + tiny));

    dF2_di = B_inos * du_inos_di - tau*beta_inos;
    dF3_da = B_arg1 * du_arg1_da - tau*beta_arg1;

    % Smooth max
    mx = max([R1h R2h],[],2);
    e1 = exp(kappa*(R1h - mx));
    e2 = exp(kappa*(R2h - mx));
    Z  = e1 + e2 + tiny;
    w1 = e1./Z;
    w2 = e2./Z;

    Dk   = mx + (1/kappa)*log(Z);
    R1h_q = nu * R1.^(nu-1) .* R1_q;
    R2h_q = nu * R2.^(nu-1) .* R2_q;
    Dk_q  = w1.*R1h_q + w2.*R2h_q;

    v_inos = R1h ./ (Dk + tiny);
    v_arg1 = R2h ./ (Dk + tiny);

    v_inos_q = (R1h_q .* Dk - R1h .* Dk_q) ./ (Dk + tiny).^2;
    v_arg1_q = (R2h_q .* Dk - R2h .* Dk_q) ./ (Dk + tiny).^2;

    % Arg_int equation diag derivative
    C4i = tau * i .* inos_max * mu_inos * q_max;
    D4i = q_max * K_inos + q_max^2 * q;

    C4a = tau * a .* arg1_max * mu_arg1 * q_max;
    D4a = q_max * K_arg1 + q_max^2 * q;

    term4i = -C4i .* ( v_inos_q .* (q ./ (D4i + tiny)) ...
                     + v_inos   .* (q_max*K_inos) ./ (D4i + tiny).^2 );

    term4a = -C4a .* ( v_arg1_q .* (q ./ (D4a + tiny)) ...
                     + v_arg1   .* (q_max*K_arg1) ./ (D4a + tiny).^2 );

    dF4_dq = term4i + term4a - tau*gamma_q;

    dF5_dn = -tau * gamma_NO * ones(N,1);
    dF6_do = -tau * gamma_orn * ones(N,1);

    div = dF1_dc + dF2_di + dF3_da + dF4_dq + dF5_dn + dF6_do;
end

function S = SourceTerm(P, x_b, sigma_hat, A_fun, t_dim, params)
    A_dim = A_fun(t_dim);
    A_hat = params.tau * A_dim;
    diff  = (P - x_b) ./ sigma_hat;
    dist2 = sum(diff.^2, 2);
    S     = A_hat * exp(-0.5*dist2);
end

function d_hat = DeathTerm(P, params, NO_max)
delta_exit = 0.02;
delta_apo  = 0.05;
K_NO       = 0.75;
n_hill     = 2;

    NO_hat = P(:,5);
    NO_dim = NO_max.*NO_hat;
    NO_pow = NO_dim.^n_hill;
    d_dim  = delta_exit + delta_apo.*(NO_pow./(K_NO^n_hill+NO_pow));
    d_hat  = params.tau .* d_dim;
end

function [P_out, w_out, V_out] = adaptive_particle_management(P_in, w_in, V_in, opts)
    [N,dim] = size(P_in);
    if dim~=6
        error('APM assumes particle states are N x 6');
    end

    w_in = w_in(:).';
    V_in = V_in(:).';

    Wtot = sum(w_in);
    if Wtot <= 0
        P_out = P_in;
        w_out = w_in(:);
        V_out = V_in(:);
        return;
    end

    N_target   = opts.N_target;
    split_fac  = opts.split_factor;
    merge_fac  = opts.merge_factor;
    d_max      = opts.d_max;
    jitter_amp = opts.split_jitter;

    Wopt = Wtot / N_target;

    % splitting
    heavy = find(w_in > split_fac*Wopt);
    P_split = [];
    w_split = [];
    V_split = [];

    for idx = heavy
        w_big = w_in(idx);
        V_big = V_in(idx);
        P0    = P_in(idx,:);

        w_child = w_big / 2;
        V_child = V_big / 2;

        jitter = jitter_amp * randn(1,dim);
        P_split = [P_split;
                   max(0,min(1,P0 + jitter));
                   max(0,min(1,P0 - jitter))];
        w_split = [w_split, w_child, w_child];
        V_split = [V_split, V_child, V_child];
    end

    survivor_mask = true(N,1);
    survivor_mask(heavy) = false;

    P_after = [P_in(survivor_mask,:); P_split];
    w_after = [w_in(survivor_mask),   w_split];
    V_after = [V_in(survivor_mask),   V_split];

    % merging
    P_tmp = P_after;
    w_tmp = w_after;
    V_tmp = V_after;

    N_now = size(P_tmp,1);
    merge_candidates = find(w_tmp < merge_fac*Wopt);

    if numel(merge_candidates) > 1
        [idxNN, distNN] = knnsearch(P_tmp(merge_candidates,:), ...
                                    P_tmp(merge_candidates,:), 'K', 2);
        neigh = idxNN(:,2);
        dist  = distNN(:,2);

        [dist_sorted, order] = sort(dist);

        used_local  = false(numel(merge_candidates),1);
        new_P_list  = [];
        new_w_list  = [];
        new_V_list  = [];

        for kk = order.'
            if dist_sorted(kk) > d_max
                break;
            end
            if used_local(kk)
                continue;
            end
            jloc = neigh(kk);
            if used_local(jloc)
                continue;
            end
            ig = merge_candidates(kk);
            jg = merge_candidates(jloc);

            wi = w_tmp(ig);   wj = w_tmp(jg);
            Vi = V_tmp(ig);   Vj = V_tmp(jg);

            w_new = wi + wj;
            V_new = Vi + Vj;
            P_new = (wi*P_tmp(ig,:) + wj*P_tmp(jg,:)) / (w_new + eps);

            new_P_list = [new_P_list; P_new];
            new_w_list = [new_w_list, w_new];
            new_V_list = [new_V_list, V_new];

            used_local(kk)   = true;
            used_local(jloc) = true;
        end

        survivor_merge = true(size(P_tmp,1),1);
        consumed_global = merge_candidates(used_local);
        survivor_merge(consumed_global) = false;

        P_out = [P_tmp(survivor_merge,:); new_P_list];
        w_out = [w_tmp(survivor_merge),   new_w_list];
        V_out = [V_tmp(survivor_merge),   new_V_list];
    else
        P_out = P_tmp;
        w_out = w_tmp;
        V_out = V_tmp;
    end

    % NOTE: no renormalization of w_out here (we keep true mass)
    w_out = w_out(:);
    V_out = V_out(:);
end

function J = jacobian_complex(F_handle, x, t, z_scalar)
    h = 1e-20;
    d = numel(x);
    J = zeros(d, d);
    for k = 1:d
        x_c = x;
        x_c(k) = x_c(k) + 1i*h;
        F_c    = F_handle(t, x_c, z_scalar);
        F_c    = drift_fix(F_c, x_c);      % enforce row shape
        J(:,k) = imag(F_c(:)) / h;
    end
end

function A_diff = diffusion_drift_internal_vec(P, w_vec, params)
    [N, dim] = size(P);
    if dim ~= 6
        error('diffusion_drift_internal_vec: P must be N x 6');
    end

    D_int   = params.D_int;
    eps_kde = params.eps_kde;
    u_floor = params.u_floor;

    w_vec = w_vec(:);
    Wtot  = sum(w_vec);
    if Wtot <= 0
        A_diff = zeros(N, dim);
        return;
    end
    w_loc = w_vec / Wtot;
    w_row = w_loc.';

    A_diff    = zeros(N, dim);
    blockSize = 200;

    for i1 = 1:blockSize:N
        i2 = min(N, i1 + blockSize - 1);
        Pi = P(i1:i2,:);
        B  = size(Pi, 1);

        DX = reshape(Pi, [B 1 dim]) - reshape(P, [1 N dim]);
        r2 = sum(DX.^2, 3);
        r  = sqrt(r2);
        q  = r / eps_kde;

        mask = (q > 0) & (q < 2);

        W    = zeros(B, N);
        dWdr = zeros(B, N);

        id1 = mask & (q <= 1);
        q1  = q(id1);
        W(id1)    = 1 - 1.5*q1.^2 + 0.75*q1.^3;
        dWdr(id1) = (-3*q1 + 2.25*q1.^2) / eps_kde;

        id2 = mask & (q > 1);
        q2  = q(id2);
        W(id2)    = 0.25*(2 - q2).^3;
        dWdr(id2) = -0.75*(2 - q2).^2 / eps_kde;

        u_block = W * w_loc;

        r_safe = r;
        r_safe(~mask) = 1;

        coeff = (dWdr ./ r_safe) .* w_row;
        coeff(~mask) = 0;

        gi = zeros(B, dim);
        for kk = 1:dim
            DXk = DX(:,:,kk);
            gi(:,kk) = sum(coeff .* DXk, 2);
        end

        u_block = max(u_block, u_floor);
        A_diff(i1:i2,:) = -D_int * (gi ./ u_block);
    end
end


