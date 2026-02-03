%% Case 1: Cybernetic model with nu=1, i.e., no ultrasensitive resource allocation

clear;clc;

% x vector = [cat inos arg1 arg_int NO orn]

% time-step size
h = 0.01;
tf = 200; % final time in hours
params.epsd = 1e-8;

% intializing vectors
t = zeros(tf/h,1);       t(1)       = 0; % time
cat = zeros(tf/h,1);     cat(1)     = 0.001; 
inos = zeros(tf/h,1);    inos(1)    = 0.001;
arg1 = zeros(tf/h,1);    arg1(1)    = 0.001;
arg_int = zeros(tf/h,1); arg_int(1) = 0.001;
NO = zeros(tf/h,1);      NO(1)      = 0.001;
orn = zeros(tf/h,1);     orn(1)     = 0.001; 
arg_ext = zeros(tf/h,1); arg_ext(1) = 1;

% cytokine signals (gaussian)
ifng_peak = 24; % in hours
il4_peak =  24*3; % in hours
ifng_spread = 10; 
il4_spread = 30;
ifng = @(x) exp(-((x-ifng_peak).^2)/(2*ifng_spread*ifng_spread));
il4  = @(x) exp(-((x-il4_peak).^2)/(2*il4_spread*il4_spread));

%% Pack all enzyme parameters into one struct

% basal synthesis rates (nM·h⁻¹) Some of them alphas are going to be zero.
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

%% Pack metabolite/cytokine parameters into the same struct

% phenotypic drift parameters (if you use mu/K form) (conc. = nM and time in hours)
params.mu.cat       = 1.5;  
params.K.cat        = 10;      
params.gamma.arg_int= 0.5;       

params.mu.inos      = 1;      
params.K.inos       = 10;      
params.gamma.NO    = 0.5;     

params.mu.arg1      = 1;   
params.K.arg1       = 10;     
params.gamma.orn    = 0.5;                    

params.K_perf       = 0.02;
params.Ab           = 100; % blood levels in micro-molar
params.gamma.arg_ext= 0.02;

params.tau          = 1/params.beta.cat;

new_final_time = tf/params.tau;
h_new = h/params.tau;

% the IFN-gamma and IL4 peaks and spreads will be shifted in Non-dimensional time. So, we need new function handles
ifng_nondim = @(x) exp(-((x-(ifng_peak/params.tau)).^2)/(2*(ifng_spread/params.tau)*(ifng_spread/params.tau)));
il4_nondim  = @(x) exp(-((x-(il4_peak/params.tau)).^2)/(2*(il4_spread/params.tau)*(il4_spread/params.tau)));
d_cat_nondim = @(cat,inos,arg1,arg_int,NO,orn,x) params.alpha.cat + ((params.kE.cat*ifng_nondim(x))/(params.KE.cat+ifng_nondim(x))) + ((params.kE.cat*il4_nondim(x))/(params.KE.cat+il4_nondim(x))) - params.beta.cat*cat;
d_inos_nondim = @(cat,inos,arg1,arg_int,NO,orn,x) params.alpha.inos + (u_inos(cat,inos,arg1,arg_int,NO,orn)*params.kE.inos*ifng_nondim(x)/(params.KE.inos+ifng_nondim(x))) - params.beta.inos*inos;
d_arg1_nondim = @(cat,inos,arg1,arg_int,NO,orn,x) params.alpha.arg1 + (u_arg1(cat,inos,arg1,arg_int,NO,orn)*params.kE.arg1*il4_nondim(x)/(params.KE.arg1+il4_nondim(x))) - params.beta.arg1*arg1;

%% -----------------------------------------------------------------------------------------------
% Let's use this space for gaussian function construction for sampling of initial conditions. The output is X0 (matrix of initial conditions)
% x vector = [cat inos arg1 arg_int NO orn]
mu=[0.01 0.01 0.01 0.01 0.01 0.01];
% Defining sigmas
alpha = 0.3; 
sigma = zeros(1,6);
for i = 1:length(mu)
    sigma(i) = alpha*min(mu(i),1-mu(i));
end
% assuming the species are correlated, we define some correlations
R = eye(6); 
shrink = 0.7;

% CAT (1) <-> Arg_int (4)
R(1, 4) = -0.7*shrink; 
R(4, 1) = -0.7*shrink; 

% iNOS (2) <-> Arg1 (3) (Competition for Arg)
R(2, 3) = -0.5*shrink;
R(3, 2) = -0.5*shrink;

% iNOS (2) <-> Arg_int (4)
R(2, 4) = -0.7*shrink;
R(4, 2) = -0.7*shrink;

% iNOS (2) <-> NO (5) (Strong Product correlation)
R(2, 5) = 0.8*shrink;
R(5, 2) = 0.8*shrink;

% Arg1 (3) <-> Orn (6) (Strong Product correlation)
R(3, 6) = 0.8*shrink;
R(6, 3) = 0.8*shrink;

D = diag(sigma);
% Building the Covariance Matrix Sigma = D * R * D
Sigma = D * R * D;

% Quick Positive Semi-Definiteness check (This is optional!)
e = eig((Sigma+Sigma')/2);
if min(e) <= 0
    error('Covariance Sigma is not positive semidefinite. Adjust R (and) sigma.');
end

% Sampling with rejection to enforce [0,1]^{6}
Np     = 1e5;         % number of particles
dim    = 6;
X0     = zeros(Np, dim);
have   = 0;

rng(1); % reproducibility

batch  = max(1000, ceil(0.02*Np));  % draw in batches. The 'ceil' function ensures that atleast 2% of Np or 1000 (whichever is smaller) samples are extracted everytime for speed
while have < Np
    draw  = mvnrnd(mu, Sigma, batch);         % candidates
    inBox = all(draw >= 0 & draw <= 1, 2);    % keep only those inside
    keep  = draw(inBox,:);
    nk    = min(size(keep,1), Np - have);
    if nk > 0
        X0(have+1:have+nk, :) = keep(1:nk, :);
        have = have + nk;
    end
end
% ---------------------------------- we get X0 --------------------------------------

P_current = max(0, min(1, X0)); % clipping values to the domain
z_current = arg_ext(1);         % initial condition on arg_ext


%% -------------------Let's define all the function handles now.---------------------
% the 'params' structure already encodes all info about the parameter values we wanna set 
% the new non-dimensional equations have some constants like: cat_max, inos_max, arg1_max, arg_int_max, NO_max, orn_max, arg_ext_max
% Let's define these constants first (see paper Appendix A).

cat_max     = max(cat(1), (params.alpha.cat+2*(params.kE.cat/(params.KE.cat+1)))/params.beta.cat);
inos_max    = max(inos(1), (params.alpha.inos+((params.kE.inos)/(params.KE.inos+1)))/params.beta.inos);
arg1_max    = max(arg1(1), (params.alpha.arg1+((params.kE.arg1)/(params.KE.arg1+1)))/params.beta.arg1);
arg_ext_max = max(arg_ext(1), (params.K_perf*params.Ab)/(params.K_perf+params.gamma.arg_ext));
arg_int_max = max(arg_int(1), (params.mu.cat*cat_max*arg_ext_max)/(params.gamma.arg_int*(params.K.cat+arg_ext_max)));
NO_max      = max(NO(1), (params.mu.inos*inos_max*arg_int_max)/(params.gamma.NO*(params.K.inos+arg_int_max)));
orn_max     = max(orn(1),(params.mu.arg1*arg1_max*arg_int_max)/(params.gamma.orn*(params.K.arg1+arg_int_max))); 

% function handles for cybernetic weights
u_inos = @(cat,inos,arg1,arg_int,NO,orn) ((inos.*inos_max.*params.mu.inos.*arg_int*arg_int_max)./(params.K.inos+arg_int*arg_int_max)) ./ max(params.epsd,((inos.*inos_max.*params.mu.inos.*arg_int*arg_int_max)./(params.K.inos+arg_int*arg_int_max) + (arg1.*arg1_max.*params.mu.arg1.*arg_int*arg_int_max)./(params.K.arg1+arg_int*arg_int_max)));
u_arg1 = @(cat,inos,arg1,arg_int,NO,orn) ((arg1.*arg1_max.*params.mu.arg1.*arg_int*arg_int_max)./(params.K.arg1+arg_int*arg_int_max)) ./ max(params.epsd,((inos.*inos_max.*params.mu.inos.*arg_int*arg_int_max)./(params.K.inos+arg_int*arg_int_max) + (arg1.*arg1_max.*params.mu.arg1.*arg_int*arg_int_max)./(params.K.arg1+arg_int*arg_int_max)));
v_inos = @(cat,inos,arg1,arg_int,NO,orn) ((inos.*inos_max.*params.mu.inos.*arg_int*arg_int_max)./(params.K.inos+arg_int*arg_int_max)) ./ max(params.epsd,max((inos.*inos_max.*params.mu.inos.*arg_int*arg_int_max)./(params.K.inos+arg_int*arg_int_max) , (arg1.*arg1_max.*params.mu.arg1.*arg_int*arg_int_max)./(params.K.arg1+arg_int*arg_int_max)));
v_arg1 = @(cat,inos,arg1,arg_int,NO,orn) ((arg1.*arg1_max.*params.mu.arg1.*arg_int*arg_int_max)./(params.K.arg1+arg_int*arg_int_max)) ./ max(params.epsd,max((inos.*inos_max.*params.mu.inos.*arg_int*arg_int_max)./(params.K.inos+arg_int*arg_int_max) , (arg1.*arg1_max.*params.mu.arg1.*arg_int*arg_int_max)./(params.K.arg1+arg_int*arg_int_max)));

% --- Vectorized Function Handles ---
% Input P is a Np X dim matrix ('dim' states, 'Np' particles)
% The output F_matrix must be Np X dim

F_particles_vec = @(t, P, z_scalar) [
% [cat,inos,arg1,arg_int,NO,orn] := [P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)]
    (params.tau*params.alpha.cat)/cat_max + (params.tau/cat_max)*(((params.kE.cat*ifng_nondim(t))/(params.KE.cat+ifng_nondim(t)))+((params.kE.cat*il4_nondim(t))/(params.KE.cat+il4_nondim(t)))) - P(:,1),...
    (params.tau*params.alpha.inos)/inos_max + (params.tau/inos_max)*((params.kE.inos*ifng_nondim(t)*u_inos(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)))/(params.KE.inos+ifng_nondim(t))) - params.beta.inos*P(:,2)*params.tau,...
    (params.tau*params.alpha.arg1)/arg1_max + (params.tau/arg1_max)*((params.kE.arg1*il4_nondim(t)*u_arg1(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)))/(params.KE.arg1+il4_nondim(t))) - params.beta.arg1*P(:,3)*params.tau,...
    (params.tau*cat_max*P(:,1)*params.mu.cat*z_scalar*arg_ext_max)./(arg_int_max*params.K.cat+arg_int_max*arg_int_max*z_scalar) - (params.tau*v_inos(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,2).*P(:,4)*inos_max*params.mu.inos*arg_int_max)./(arg_int_max*params.K.inos+arg_int_max*arg_int_max*P(:,4)) - (params.tau*v_arg1(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,3).*P(:,4)*arg1_max*params.mu.arg1*arg_int_max)./(arg_int_max*params.K.arg1+arg_int_max*arg_int_max*P(:,4)) - params.tau*params.gamma.arg_int*P(:,4),...
    (params.tau*v_inos(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,2).*P(:,4)*inos_max*params.mu.inos*arg_int_max)./(NO_max*params.K.inos+NO_max*arg_int_max*P(:,4)) - params.tau*params.gamma.NO*P(:,5),...
    (params.tau*v_arg1(P(:,1),P(:,2),P(:,3),P(:,4),P(:,5),P(:,6)).*P(:,3).*P(:,4)*arg1_max*params.mu.arg1*arg_int_max)./(orn_max*params.K.arg1+orn_max*arg_int_max*P(:,4)) - params.tau*params.gamma.orn*P(:,6)
  ]; 
   
   
% --- initial weights -------
w_vec = (1/Np) * ones(Np, 1);
%----------------------------


% ---- source term inputs-------------------------------------------------------------------------
x_b = [0 0 0 0 0 0];        % baseline
sigma_hat = 0.04*ones(1,6); % widths 
A_fun = @(t) (0.02 + 0.01*sin(0.1*t)); % this can be changed to match the monocyte influx rates
% ------------------------------------------------------------------------------------------------


% -------- parameters for APM --------------------------------------------------------------------
APM_opts.N_target     = Np;       % keep resolution equal to initial
APM_opts.merge_factor = 0.5;      % triggers merging below 0.5*Wopt
APM_opts.split_factor = 2.0;      % triggers splitting above 2*Wopt
APM_opts.d_max        = 0.02;     % merging distance tolerance
APM_opts.split_jitter = 0.005;    % small jitter for KDE support

% We can tune:
% d_max too small → little merging
% d_max too large → over-smoothing
% jitter too large → artificial diffusion
% split_factor lower → more splitting
% ------------------------------------------------------------------------------------------------

% ------------------------- parameters for 'aux' struct in analytic divergence -------------------
aux.inos_max    = inos_max;
aux.arg1_max    = arg1_max;
aux.arg_int_max = arg_int_max;
aux.ifng_nondim = ifng_nondim;
aux.il4_nondim  = il4_nondim;
aux.kappa       = 20;   % or 10, 30, tune as needed
% ------------------------------------------------------------------------------------------------


% external arginine evolution (scalar field) function handle 
R_ext_handle       = @(z_scalar) ((params.K_perf*params.tau)/arg_ext_max)*(params.Ab - z_scalar*arg_ext_max) - params.gamma.arg_ext*z_scalar*params.tau;
uptake_rate_handle = @(P, z_scalar) (params.tau*P(:,1)*cat_max*params.mu.cat*arg_ext_max*z_scalar)./(arg_int_max*params.K.cat+arg_int_max*arg_ext_max*z_scalar);
% The External Metabolite Derivative Function: dz/dt
F_external_z = @(t, P, z_scalar, w_vec) -sum(w_vec(:) .* uptake_rate_handle(P, z_scalar)) + R_ext_handle(z_scalar);


%% ------------------------ Vectorized SSP-RK3 (Time Loop) ------------------------
% Setup: P_current, z_current, dt, w_vec, etc. initialized as before
tn = t(1);
N_steps = round(tf / h); 
dt = h_new;
t_save = [4, 20, 40, 80, 120, 160, 195]/params.tau; % non-dimensional times for KDE snapshots
Nt = numel(t_save);

inos_x = cell(1,Nt); % --
inos_f = cell(1,Nt); %  | storage for KDE grids and PDFs
arg1_x = cell(1,Nt); %  |
arg1_x = cell(1,Nt); % --

cloud_inos_arg1 = cell(1,Nt); % storage for cluster visualization
cloud_w         = cell(1,Nt); % "

inos_mean = nan(1,Nt); % for storing mean values of the marginal distribution
arg1_mean = nan(1,Nt);

%profile on
for n = 1:N_steps
    % ----------------------------------------------------
    % STAGE 1: Derivative at t_n
    % ----------------------------------------------------
    kP_n = F_particles_vec(tn, P_current, z_current);
	kw_n = -(w_vec(:).*analytic_divergence(tn, P_current, z_current, params, aux)) + SourceTerm(P_current, x_b, sigma_hat, A_fun, tn, params); %w_vec(:) prevents broadcasting
    kz_n = F_external_z(tn, P_current, z_current, w_vec);
    
    % Evolve to intermediate state 1
    P1 = P_current + dt * kP_n;
	w1 = w_vec + dt * kw_n; w1=w1(:);
    z1 = z_current + dt * kz_n;

    % ----------------------------------------------------
    % STAGE 2: Derivative at t_{n+1} (tn_plus_dt)
    % ----------------------------------------------------
    tn_plus_dt = tn + dt;
    % Calculate derivatives based on intermediate state 1 (P1, z1)
    kP_1 = F_particles_vec(tn_plus_dt, P1, z1);
	kw_1 = -(w1(:).*analytic_divergence(tn_plus_dt, P1, z1, params, aux)) + SourceTerm(P1, x_b, sigma_hat, A_fun, tn_plus_dt, params); 
    kz_1 = F_external_z(tn_plus_dt, P1, z1, w1);
    
    % Evolve to intermediate state 2 (SSP averaging step)
    P2 = (3/4)*P_current + (1/4)*P1 + (1/4)*dt * kP_1;
	w2 = (3/4)*w_vec + (1/4)*w1 + (1/4)*dt * kw_1;w2=w2(:);
    z2 = (3/4)*z_current + (1/4)*z1 + (1/4)*dt * kz_1;
    
    % ----------------------------------------------------
    % STAGE 3 (Final): Derivative at t_{n+1} (tn_plus_dt)
    % ----------------------------------------------------
    % Calculate derivatives based on intermediate state 2 (P2, z2)
    % *** We correctly reuse the time point t_n + dt ***
    kP_2 = F_particles_vec(tn_plus_dt, P2, z2); 
	kw_2 = -(w2(:).*analytic_divergence(tn_plus_dt, P2, z2, params, aux)) + SourceTerm(P2, x_b, sigma_hat, A_fun, tn_plus_dt, params); 
    kz_2 = F_external_z(tn_plus_dt, P2, z2, w2);
    
    % Final SSP averaging step for t(n+1)
    P_next = (1/3)*P_current + (2/3)*P2 + (2/3)*dt * kP_2;
	w_next = (1/3)*w_vec + (2/3)*w2 + (2/3)*dt * kw_2;w_next = w_next(:);
    z_next = (1/3)*z_current + (2/3)*z2 + (2/3)*dt * kz_2;

    % ----------------------------------------------------
    % Update and Housekeeping
    % ----------------------------------------------------
    P_current = max(0, P_next); 
    z_current = max(0, z_next); 
	w_vec     = max(0, w_next(:));

 % -------------------------------------------------------
    w_vec = max(0, w_next(:));                              %|
Wtot  = sum(w_vec);                                         %|
if Wtot <= 0                                                %|
    error('All weights vanished.');                         %| this block is to help with the blowup in weights.
end                                                         %|
w_vec = w_vec / Wtot;  % probability measure, sum(w) = 1    %|
 % -------------------------------------------------------

    tn = tn_plus_dt; % tn is advanced by dt
	
	% ------------------KDE snapshots 1D --------------------
	w_norm = w_vec(:)/sum(w_vec); % normalize weights
	tol_t = 0.5*dt; % tolerance so tn 'matches' a requested save time
	idx_hit = find(abs(tn - t_save) <= tol_t, 1, 'first');
	if ~isempty(idx_hit)
	    k = idx_hit;
		inos_vals = P_current(:,2); % pull inos coordinate 
		arg1_vals = P_current(:,3); % pull arg1 coordinate

        fprintf('t = %.2f (nd), std(iNOS) = %.3e, min = %.3e, max = %.3e\n', ...  % diagnostics 
        tn, std(inos_vals), min(inos_vals), max(inos_vals));                      %  "
        fprintf('   sum(w) = %.6e, min(w) = %.3e, max(w) = %.3e\n', ...           %  "
        sum(w_vec), min(w_vec), max(w_vec));                                      %  "

		[f_inos, x_inos] = ksdensity(inos_vals, 'Weights', w_norm,'NumPoints',200,'Support',[0 1],'BoundaryCorrection','reflection', 'Bandwidth',0.02); % weighted KDE
		[f_arg1, x_arg1] = ksdensity(arg1_vals, 'Weights', w_norm,'NumPoints',200,'Support',[0 1],'BoundaryCorrection','reflection', 'Bandwidth',0.02);

        % ----- store particle cloud for visualization (subsample) -----
        Nsamp = min(2000, size(P_current,1));
        idx_samp = randperm(size(P_current,1), Nsamp);
        cloud_inos_arg1{k} = P_current(idx_samp, [2 3]);  % [iNOS, Arg1]
        cloud_w{k}         = w_norm(idx_samp);

        fprintf('   integral pdf ~ %.3e\n', trapz(x_inos, f_inos)); % --> diagnostics

		inos_x{k}      = x_inos;
        inos_f{k}      = f_inos;
        arg1_x{k}      = x_arg1;
        arg1_f{k}      = f_arg1;
        dx_inos        = x_inos(2) - x_inos(1);
        iNOS_mean(k)   = sum(x_inos(:) .* f_inos(:)) * dx_inos;
        dx_arg1        = x_arg1(2) - x_arg1(1);
        arg1_mean(k)   = sum(x_arg1(:) .* f_arg1(:)) * dx_arg1;
        % ---------- NEW: store particle cloud in (iNOS, Arg1) ----------
    Nsamp    = min(2000, size(P_current,1));   % subsample for plotting
    idx_samp = randperm(size(P_current,1), Nsamp);

    cloud_inos_arg1{k} = P_current(idx_samp, [2 3]);  % [iNOS, Arg1]
    cloud_w{k}         = w_norm(idx_samp);            % normalized weights
	end
	
	% ------------------- APM ------------------
	% if mod(n,20) == 0
   %  [P_current, w_vec] = adaptive_particle_management(P_current, w_vec, APM_opts);
    % end
end
%profile viewer



%% --------------------------------------------- Results and Visualizations -------------------------------------------------
% 1)  e.g. show iNOS distribution snapshots
figure;
for k = 1:Nt
    if isempty(inos_x{k}), continue; end
    plot(inos_x{k}, inos_f{k}, 'DisplayName', sprintf('t = %.1f h', t_save(k)));
    hold on;
end
xlabel('iNOS (nondimensional)');
ylabel('PDF');
title('iNOS KDE snapshots');
legend show;

figure;
for k = 1:Nt
    if isempty(arg1_x{k}), continue; end
    plot(arg1_x{k}, arg1_f{k}, 'DisplayName', sprintf('t = %.1f h', t_save(k)));
    hold on;
end
xlabel('arg1 (nondimensional)');
ylabel('PDF');
title('arg1 KDE snapshots');
legend show;

% 2) Visualizing actual clustering in state space
figure;
for k = 1:Nt
    if isempty(cloud_inos_arg1{k}), continue; end

    subplot(2, ceil(Nt/2), k);
    pts = cloud_inos_arg1{k};
    scatter(pts(:,1), pts(:,2), 8, cloud_w{k}, 'filled');  % color = weight
    xlabel('iNOS (nd)');
    ylabel('Arg1 (nd)');
    title(sprintf('t = %.1f (non-dim)', t_save(k)));
    xlim([0 1]); ylim([0 1]);
    colorbar;
end
sgtitle('Particle cloud in (iNOS, Arg1) space');

%============== improved visualization =====================
idx_subset = [1 2 4 7];   % t = 0.8, 4.0, 16.0, 39.0 (non-dim)

figure('Units','centimeters','Position',[2 2 14 10]);
tl = tiledlayout(2,2,'TileSpacing','none','Padding','compact');

for j = 1:numel(idx_subset)
    k = idx_subset(j);
    if isempty(cloud_inos_arg1{k}), continue; end

    ax = nexttile(tl);

    pts   = cloud_inos_arg1{k};   % [iNOS, Arg1]
    cvals = cloud_w{k};           % weights

    scatter(ax, pts(:,1), pts(:,2), 20, cvals, ...
            'filled','MarkerEdgeColor','none');
    hold(ax,'on');

    % one colorbar per tile
    cb = colorbar(ax);
    cb.Label.String = 'Particle weight';
    cb.FontSize     = 7;

    % ---- auto-zoom around the cloud ----
    x = pts(:,1);  y = pts(:,2);
    xr = [min(x) max(x)];
    yr = [min(y) max(y)];
    if xr(1) == xr(2), xr = xr + [-1 1]*1e-3; end
    if yr(1) == yr(2), yr = yr + [-1 1]*1e-3; end

    xpad = 0.2 * max(1e-3, xr(2) - xr(1));
    ypad = 0.2 * max(1e-3, yr(2) - yr(1));

    xmin = max(0, xr(1) - xpad);
    xmax = min(1, xr(2) + xpad);
    ymin = max(0, yr(1) - ypad);
    ymax = min(1, yr(2) + ypad);

    xlim(ax,[xmin xmax]);
    ylim(ax,[ymin ymax]);

    axis(ax,'square');
    set(ax,'FontSize',8,'LineWidth',1);

    title(ax, sprintf('t = %.1f (non-dim)', t_save(k)), ...
          'FontWeight','normal');

    if j > 2
        xlabel(ax,'iNOS (nd)');
    end
    if mod(j-1,2) == 0
        ylabel(ax,'Arg1 (nd)');
    end
end

sgtitle(tl,'Contraction of particle cloud in (iNOS, Arg1) space', ...
        'FontWeight','bold','FontSize',9);
%==============================================================================   


% 3) Overlay the vector field (phase portrait) at a representative time
% by Choosing a representative time and external arginine
t_phase = t_save(end);       % last snapshot time (nd)
z_phase = z_current;         % or store z(t) in an array during the run

% Fix other coordinates at something reasonable (e.g. final means)
cat0    = mean(P_current(:,1));
q0      = mean(P_current(:,4));   % Arg_int
NO0     = mean(P_current(:,5));
orn0    = mean(P_current(:,6));

% Build a grid in (iNOS, Arg1)
[inos_grid, arg1_grid] = meshgrid(linspace(0,1,20), linspace(0,1,20));
Ng = numel(inos_grid);

P_grid = zeros(Ng, 6);
P_grid(:,1) = cat0;
P_grid(:,2) = inos_grid(:);
P_grid(:,3) = arg1_grid(:);
P_grid(:,4) = q0;
P_grid(:,5) = NO0;
P_grid(:,6) = orn0;

F_grid = F_particles_vec(t_phase, P_grid, z_phase);

u = reshape(F_grid(:,2), size(inos_grid));  % d(iNOS)/dt
v = reshape(F_grid(:,3), size(arg1_grid));  % d(Arg1)/dt

figure;
quiver(inos_grid, arg1_grid, u, v); hold on;

% Overlay final particle cloud
pts_final = cloud_inos_arg1{end};
scatter(pts_final(:,1), pts_final(:,2), 10, 'r', 'filled');

xlabel('iNOS (nd)');
ylabel('Arg1 (nd)');
title(sprintf('Phase portrait in (iNOS, Arg1) at t = %.1f h', t_phase));
xlim([0 1]); ylim([0 1]);
legend('vector field','particles');

% 4)  COMPUTE ATTRACTOR AT FINAL SNAPSHOT TIME
w_norm_final = w_vec(:) / sum(w_vec);     % normalize weights
x_star = sum(P_current .* w_norm_final, 1);   % 1×6 row vector
z_star = z_current;                         % final external arginine level
t_star = t_save(end);                       % the time of the last snapshot

fprintf('\n===== ATTRACTOR CANDIDATE x_star =====\n');
disp(x_star);

%% --- Evaluate F(x_star) to check if it is an equilibrium ---
Fx = F_particles_vec(t_star, x_star, z_star);   % 1×6
Fx_norm = norm(Fx);

fprintf('||F(x_star)|| = %.3e\n', Fx_norm);
if Fx_norm < 1e-6
    disp('=> x_star is approximately a FIXED POINT of the 6-D system.');
else
    disp('=> x_star is NOT a fixed point; dynamics approach a trajectory.');
end

%% Compute Jacobian
J = jacobian_complex(F_particles_vec, x_star, t_star, z_star);

fprintf('\n===== EIGENVALUES OF JACOBIAN @ ATTRACTOR =====\n');
eigvals = eig(J);
disp(eigvals);

% Check stability
if all(real(eigvals) < 0)
    disp('=> The attractor is a LOCALLY STABLE sink.');
elseif any(real(eigvals) > 0)
    disp('=> The attractor is UNSTABLE or a saddle.');
else
    disp('=> The attractor has neutrally stable directions.');
end
% ----------------------------------------------------------------------------------------------------------------------------


% complex -step differentiation is slowing down the performance due to 18 vectorized calls. so we calculated the exact analytic divergence by replacing max(a,b) with 'logsumexp'
% here is the function for the analytic divergence
function div = analytic_divergence(t, P, z_scalar, params, aux)
% ANALYTIC_DIVERGENCE
%   Compute divergence ∇·F(t,x,z) for each particle x in P.
%
% INPUTS:
%   t        : scalar time (nondimensional)
%   P        : N x 6 particle states [cat, iNOS, Arg1, Arg_int, NO, orn]
%   z_scalar : scalar external arginine (nondimensional) (not used in divergence)
%   params   : struct with fields like params.tau, params.mu.*, params.K.*,
%              params.beta.*, params.gamma.*, params.kE.*, params.KE.*
%   aux      : struct with:
%                aux.inos_max      (inos_max)
%                aux.arg1_max      (arg1_max)
%                aux.arg_int_max   (arg_int_max)
%                aux.ifng_nondim   (handle: @(t)...)
%                aux.il4_nondim    (handle: @(t)...)
%                aux.kappa         (smooth-max sharpness)
%
% OUTPUT:
%   div      : N x 1 vector, divergence at each particle

    % unpack sizes and state
    [N,dim] = size(P);
    if dim ~= 6
        error('analytic_divergence: P must be N x 6.');
    end

    c = P(:,1);  % cat
    i = P(:,2);  % iNOS
    a = P(:,3);  % Arg1
    q = P(:,4);  % Arg_int
    n = P(:,5);  %#ok<NASGU> % NO  (unused in derivatives)
    o = P(:,6);  %#ok<NASGU> % Orn (unused in derivatives)

    % ----- unpack parameters -----
    tau   = params.tau;

    % enzyme "max" scales
    inos_max    = aux.inos_max;
    arg1_max    = aux.arg1_max;
    q_max       = aux.arg_int_max;

    % rate parameters
    mu_inos     = params.mu.inos;
    mu_arg1     = params.mu.arg1;

    K_inos      = params.K.inos;
    K_arg1      = params.K.arg1;

    beta_inos   = params.beta.inos;
    beta_arg1   = params.beta.arg1;

    gamma_q     = params.gamma.arg_int;  % Arg_int
    gamma_NO    = params.gamma.NO;
    gamma_orn   = params.gamma.orn;

    % kinetic params for inducible synthesis
    kE_inos     = params.kE.inos;
    kE_arg1     = params.kE.arg1;

    KE_inos     = params.KE.inos;
    KE_arg1     = params.KE.arg1;

    % smooth-max sharpness
    kappa       = aux.kappa;

    % cytokine signals (already nondimensional)
    ifng_val    = aux.ifng_nondim(t);
    il4_val     = aux.il4_nondim(t);

    % tiny regularizer to avoid 0/0
    tiny = 1e-12;

    % =====================================
    % 1) DIAGONAL DERIVATIVE FOR F1: cat
    % -------------------------------------
    % F1 = constant(t) - c  → ∂F1/∂c = -1
    dF1_dc = -ones(N,1);

    % ===============================================================
    % 2) BUILD R1, R2 AND THEIR q- and i/a-DERIVATIVES (for u and v)
    % ---------------------------------------------------------------
    % Raw rates:
    %   R1(i,q) = (i * inos_max * mu_inos * q * q_max) / (K_inos + q*q_max)
    %   R2(a,q) = (a * arg1_max * mu_arg1 * q * q_max) / (K_arg1 + q*q_max)

    A1   = i .* inos_max * mu_inos * q_max;   % N x 1
    den1 = K_inos + q * q_max;                % N x 1
    R1   = (A1 .* q) ./ (den1 + tiny);        % N x 1

    A2   = a .* arg1_max * mu_arg1 * q_max;   % N x 1
    den2 = K_arg1 + q * q_max;                % N x 1
    R2   = (A2 .* q) ./ (den2 + tiny);        % N x 1

    % ∂R1/∂q and ∂R2/∂q
    R1_q = A1 .* K_inos  ./ (den1 + tiny).^2;
    R2_q = A2 .* K_arg1  ./ (den2 + tiny).^2;

    % For u-weights we need C1(q), C2(q) = ∂R1/∂i, ∂R2/∂a
    %   C1(q) = dR1/di = inos_max * mu_inos * q*q_max / (K_inos + q*q_max)
    %   C2(q) = dR2/da = arg1_max * mu_arg1 * q*q_max / (K_arg1 + q*q_max)
    C1 = (inos_max * mu_inos * q .* q_max) ./ (den1 + tiny);   % N x 1
    C2 = (arg1_max * mu_arg1 * q .* q_max) ./ (den2 + tiny);   % N x 1

    % Sum S = R1 + R2 for u-weights
    S  = R1 + R2 + tiny;

    % =========================================
    % 3) u-WEIGHTS AND THEIR i/a-DERIVATIVES
    % -----------------------------------------
    % u_inos   = R1 / (R1+R2)
    % u_arg1   = R2 / (R1+R2)
    %
    % du_inos/di   = C1 * R2 / S^2
    % du_arg1/da   = C2 * R1 / S^2

    u_inos       = R1 ./ S;
    u_arg1       = R2 ./ S;

    du_inos_di   = (C1 .* R2) ./ (S.^2);
    du_arg1_da   = (C2 .* R1) ./ (S.^2);

    % ==================================================
    % 4) DIAGONAL DERIVATIVE FOR F2 (iNOS) AND F3 (Arg1)
    % --------------------------------------------------
    % F2 = τ α_inos / inos_max + B_inos(t) * u_inos - τ β_inos * i
    %      where B_inos(t) = τ/inos_max * kE_inos * IFN(t) / (KE_inos + IFN(t))
    %
    % dF2/di = B_inos(t) * du_inos/di - τ β_inos
    %
    % F3 = τ α_arg1 / arg1_max + B_arg1(t) * u_arg1 - τ β_arg1 * a
    %      where B_arg1(t) = τ/arg1_max * kE_arg1 * IL4(t) / (KE_arg1 + IL4(t))
    %
    % dF3/da = B_arg1(t) * du_arg1/da - τ β_arg1

    B_inos  = tau / inos_max * (kE_inos * ifng_val / (KE_inos + ifng_val + tiny));
    B_arg1  = tau / arg1_max * (kE_arg1 * il4_val / (KE_arg1 + il4_val + tiny));

    dF2_di  = B_inos * du_inos_di - tau * beta_inos;
    dF3_da  = B_arg1 * du_arg1_da - tau * beta_arg1;

    % ==============================================================
    % 5) SMOOTH-MAX AND v-WEIGHTS (for F4, F5, F6 derivatives wrt q)
    % --------------------------------------------------------------
    % Smooth max D_kappa = (1/kappa)*log(exp(kappa R1) + exp(kappa R2))
    % Use stable version: subtract max before exponentiating.
    %
    % w1 = softmax weight for R1, w2 = for R2
    %
    % D_kappa_q = w1 R1_q + w2 R2_q
    %
    % v_inos   = R1 / D_kappa
    % v_arg1   = R2 / D_kappa
    %
    % v_inos_q = (R1_q D_kappa - R1 D_kappa_q) / D_kappa^2
    % v_arg1_q = (R2_q D_kappa - R2 D_kappa_q) / D_kappa^2

    mx   = max([R1 R2], [], 2);          % N x 1
    e1   = exp(kappa * (R1 - mx));
    e2   = exp(kappa * (R2 - mx));
    Z    = e1 + e2 + tiny;

    w1   = e1 ./ Z;
    w2   = e2 ./ Z;

    Dk   = mx + (1/kappa) * log(Z);      % smooth max
    Dk_q = w1 .* R1_q + w2 .* R2_q;

    v_inos   = R1 ./ (Dk + tiny);
    v_arg1   = R2 ./ (Dk + tiny);

    v_inos_q = (R1_q .* Dk - R1 .* Dk_q) ./ (Dk + tiny).^2;
    v_arg1_q = (R2_q .* Dk - R2 .* Dk_q) ./ (Dk + tiny).^2;

    % ===========================================
    % 6) DIAGONAL DERIVATIVE FOR F4: Arg_int (q)
    % -------------------------------------------
    % F4 = (cat-intake term, independent of q)
    %    - C4i * v_inos * q / D4i(q)
    %    - C4a * v_arg1 * q / D4a(q)
    %    - τ γ_arg_int q
    %
    % where:
    %   C4i = τ * i * inos_max * μ_inos * q_max
    %   D4i(q) = q_max*K_inos + q_max^2*q
    %
    %   C4a = τ * a * arg1_max * μ_arg1 * q_max
    %   D4a(q) = q_max*K_arg1 + q_max^2*q
    %
    % For G(q) = C * v(q) * q / D(q) with D(q)=A+Bq:
    %   dG/dq = C [ v_q * q / D + v * A / D^2 ]
    %
    % So:
    %   dF4/dq = -C4i [ v_inos_q * q / D4i + v_inos * (q_max*K_inos)/D4i^2 ]
    %            -C4a [ v_arg1_q * q / D4a + v_arg1 * (q_max*K_arg1)/D4a^2 ]
    %            - τ γ_arg_int

    C4i = tau * i .* inos_max * mu_inos * q_max;
    D4i = q_max * K_inos  + q_max^2 * q;

    C4a = tau * a .* arg1_max * mu_arg1 * q_max;
    D4a = q_max * K_arg1  + q_max^2 * q;

    term4i = -C4i .* ( v_inos_q .* (q ./ (D4i + tiny)) ...
                     + v_inos   .* (q_max * K_inos) ./ (D4i + tiny).^2 );

    term4a = -C4a .* ( v_arg1_q .* (q ./ (D4a + tiny)) ...
                     + v_arg1   .* (q_max * K_arg1) ./ (D4a + tiny).^2 );

    dF4_dq = term4i + term4a - tau * gamma_q;

    % =================================================
    % 7) DIAGONAL DERIVATIVES FOR F5 (NO) AND F6 (Orn)
    % -------------------------------------------------
    % F5 = [production term v_inos(i,a,q)] - τ γ_NO * n
    %   → ∂F5/∂n = - τ γ_NO
    %
    % F6 = [production term v_arg1(i,a,q)] - τ γ_orn * o
    %   → ∂F6/∂o = - τ γ_orn

    dF5_dn = -tau * gamma_NO;
    dF6_do = -tau * gamma_orn;

    dF5_dn = dF5_dn * ones(N,1);
    dF6_do = dF6_do * ones(N,1);

    % =========================================================================================
    % 8) FINAL DIVERGENCE ∇·F
    % -----------------------------------------------------------------------------------------
    % ∇·F = ∂F1/∂c + ∂F2/∂i + ∂F3/∂a + ∂F4/∂q + ∂F5/∂n + ∂F6/∂o

    div = dF1_dc + dF2_di + dF3_da + dF4_dq + dF5_dn + dF6_do;
end  


%% --------------------------------- Source Term ------------------------------
function S = SourceTerm(P, x_b, sigma_hat, A_fun, t, params)
% GaussianSourceTime computes nondimensional Gaussian source for particles
%
% P         : Np x d  matrix of nondimensional particle states
% x_b       : 1 x d   nondimensional baseline state
% sigma_hat : 1 x d   nondimensional Gaussian widths
% A_fun     : function handle giving *dimensional* amplitude, A(t)
% t         : scalar time (dimensional time input to A_fun)
%
% RETURNS:
% S : Np x 1 vector of nondimensional source values

A_dim = A_fun(t); % evaluate amplitude in dimensional units
A_hat = params.tau*A_dim ; % convert amplitude to non-dim form
diff = (P - x_b)./ sigma_hat; % Np x d 
dist2 = sum(diff.^2,2); % Np x 1
S = A_hat * exp(-0.5*dist2); % evaluate gaussian
end
% ----------------------------------------------------------------------------------------------------------------------------

%% ------------------------------------------ Adaptive particle management ---------------------------------------------------
function [P_out, w_out] = adaptive_particle_management(P_in, w_in, opts)
% Inputs:
%   P_in : N x 6   particle states   (in [0,1])
%   w_in : 1 x N   particle weights
% Output:
%   P_out, w_out : updated particles & weights

% setup 
[N,dim] = size(P_in);
if dim~=6
   error('APM assumes particle states are N X 6')
end

Wtot = sum(w_in);
if Wtot <= 0
P_out = P_in;
w_out = w_in;
return;
end

N_target    = opts.N_target;
split_fac   = opts.split_factor;
merge_fac   = opts.merge_factor;
d_max       = opts.d_max;
jitter_amp  = opts.split_jitter;

Wopt = Wtot / N_target; % optimal weight scaling

% ------- splitting for overweight particles ---------
heavy = find(w_in > split_fac*Wopt);
P_split = [];
w_split = [];

for idx = heavy
  w_big   = w_in(idx);
  w_child = w_big / 2;
  P0      = P_in(idx,:);

  jitter = jitter_amp * randn(1,dim); % symmetric jitter for kernel support

  P_split = [P_split;
             max(0,min(1,P0+jitter));
             max(0,min(1,P0-jitter))];
  w_split = [w_split, w_child, w_child];
end
  
survivor_mask = true(N,1);
survivor_mask(heavy) = false;
P_after_split = [P_in(survivor_mask,:); P_split];
w_after_split = [w_in(survivor_mask)  , w_split];

% --------------- merging for tiny particles ---------------
P_tmp = P_after_split;
    w_tmp = w_after_split;

    N_now = size(P_tmp,1);
    merge_candidates = find(w_tmp < merge_fac * Wopt);

    if numel(merge_candidates) > 1

        % kd-tree like NN query (requires knnsearch)
        [idxNN, distNN] = knnsearch(P_tmp(merge_candidates,:), ...
                                   P_tmp(merge_candidates,:), 'K', 2);

        neigh = idxNN(:,2);
        dist  = distNN(:,2);

        [dist_sorted, order] = sort(dist);

        used_local  = false(numel(merge_candidates),1);
        new_P_list  = [];
        new_w_list  = [];

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

            wi = w_tmp(ig);
            wj = w_tmp(jg);

            w_new = wi + wj;
            P_new = (wi*P_tmp(ig,:) + wj*P_tmp(jg,:)) / w_new;   % p-scheme

            new_P_list = [new_P_list; P_new];
            new_w_list = [new_w_list, w_new];

            used_local(kk)  = true;
            used_local(jloc) = true;
        end

        % Remove merged-away particles
        survivor_merge = true(size(P_tmp,1),1);
        consumed_global = merge_candidates(used_local);
        survivor_merge(consumed_global) = false;

        P_out = [P_tmp(survivor_merge,:); new_P_list];
        w_out = [w_tmp(survivor_merge),   new_w_list];
    else
        % No merging possible
        P_out = P_tmp;
        w_out = w_tmp;
    end

    % ------------------------------
    % 3) Mass renormalization
    % ------------------------------
    Wnew = sum(w_out);
    if Wnew > 0
       % w_out = w_out * (Wtot / Wnew);
        w_out = w_out / Wnew; 
    end

end
% ----------------------------------------------------------------------------------------------------------------------------
%% =====================================================
%        JACOBIAN OF F AT THE ATTRACTOR (complex-step)
% ======================================================
function J = jacobian_complex(F_handle, x, t, z_scalar)
    % F_handle: function handle  (t, X, z) → F
    % x: 1×6 evaluation point
    % t, z_scalar: scalar inputs
    % Returns: 6×6 Jacobian
    h = 1e-20;                       % complex-step magnitude
    d = numel(x);
    J = zeros(d, d);
    for k = 1:d
        x_c = x;
        x_c(k) = x_c(k) + 1i*h;      % perturb kth dimension
        F_c = F_handle(t, x_c, z_scalar);
        J(:,k) = imag(F_c(:)) / h;   % exact derivative
    end
end

