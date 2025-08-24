clear
clc
rng(1);    % for reprodusibility

%% ============================ Section 2 =============================
%% 2.1 Quarterly Data: New Moment Conditions
%% ===================================================================

%% 2.1.1 a MATLAB function that returns as first output argument the given vector & matrix --> Q_cbm_moments.m

%% 2.1.2 Construct weighting matrix W_T 
N = 25;
tau = 1500;
W_matrix = [eye(N), zeros(N,1); zeros(1,N), tau];


%% ========================================================================
%% 2.2 A Full Resurrection of the CBM
%% ===================================================================

%% 2.2.1 Setup GMM and estimation inputs
gmmopt = struct();
gmmopt.infoz.momt = 'Q_cbm_moments';   % moment function file
gmmopt.W0 = 'Win';                    % custom weighting matrix
gmmopt.gmmit = 1;                     % one-step GMM
gmmopt.hess = 'gn';                   % Gauss-Newton
gmmopt.S = 'NW';                      % Newey-West spectral density
gmmopt.plot = 0;                      % no plots
gmmopt.prt = 1;                       % print results

% Initial parameter guess: [alpha, mu, gamma]
init_q = [0; 1; 3];

% Extract data and construct matrices
Q_raw = readtable('QuarterlyData_1949Q1_2024Q1.csv');
r_excess_q = table2array(Q_raw(:, 5:29)) - table2array(Q_raw(:, 4));  % realized - risk-free

% ---------- FILTERED consumption ----------
X_q_filtered = [table2array(Q_raw(:,2)), r_excess_q];
Y_q = zeros(280, 26);           % all moments - 0 
Z_q = ones(280,1);              % set to vector of ones since no instruments used

% Call Gmm for filtered data 
[est_q_filt, ~] = gmm(init_q, gmmopt, Y_q, X_q_filtered, Z_q, W_matrix);


% Package results into table
results_q_filt = table(est_q_filt.b, est_q_filt.se, ...
    est_q_filt.b - 1.96 * est_q_filt.se, est_q_filt.b + 1.96 * est_q_filt.se, ...
    'VariableNames', {'Estimate', 'StdErr', 'CI_Lower', 'CI_Upper'});

% disp(results_q_filt);   

%% 2.2.2 Hypothesis Tests for gamma ≤ 10 and alpha = 0
tstat_gamma_filt = (est_q_filt.b(3) - 10) / est_q_filt.se(3);
tstat_alpha_filt = (est_q_filt.b(1) - 0) / est_q_filt.se(1);
pval_gamma_filt = 1 - normcdf(tstat_gamma_filt);
pval_alpha_filt = 2 * (1 - normcdf(abs(tstat_alpha_filt)));

% Result Interpretation --> Please refer pdf

%% 2.2.3 Compare model-implied returns to mean observed returns
beta_q = 0.95;
sdf_q_filt = beta_q * X_q_filtered(:,1) .^ (-est_q_filt.b(3));
ret_model_q_filt = est_q_filt.b(1) - mean((sdf_q_filt - est_q_filt.b(2)) .* r_excess_q, 1) ./ est_q_filt.b(2);
ret_sample_q = mean(r_excess_q, 1);

% convert to percent
ret_model_pct_q_filt = ret_model_q_filt * 100;
ret_sample_pct_q = ret_sample_q * 100;

% Plot --> scatter for Filtered Consumption
figure;
scatter(ret_model_pct_q_filt, ret_sample_pct_q, 30, [1 0 1], 'filled');  % magenta dots
hold on;
diag_min = min([ret_model_pct_q_filt, ret_sample_pct_q]);
diag_max = max([ret_model_pct_q_filt, ret_sample_pct_q]);
plot([diag_min diag_max], [diag_min diag_max], 'k-', 'LineWidth', 0.5);  % solid thin black line
xlabel('Model-Implied Excess Return (%)');
ylabel('Realized Mean Excess Return (%)');
title('Filtered Consumption: Model vs. Data');
xlim([diag_min diag_max]); ylim([diag_min diag_max]);
box on;
hold off;

% Interpretation of plot : Please refer pdf

%% 2.2.4 Cross-sectional R square
r2_q_filt = 1 - var(ret_sample_q - ret_model_q_filt) / var(ret_sample_q);

% Interpretation --> please refer pdf


%% Repeated process --> question 2-4 for UNFILTERED consumption
X_q_unfiltered = [table2array(Q_raw(:,3)), r_excess_q];

% Run GMM
[est_q_unfilt, ~] = gmm(init_q, gmmopt, Y_q, X_q_unfiltered, Z_q, W_matrix);


% Results
results_q_unfilt = table(est_q_unfilt.b, est_q_unfilt.se, ...
    est_q_unfilt.b - 1.96 * est_q_unfilt.se, est_q_unfilt.b + 1.96 * est_q_unfilt.se, ...
    'VariableNames', {'Estimate', 'StdErr', 'CI_Lower', 'CI_Upper'});
% disp(results_q_unfilt);

% Hypothesis Tests --> the same way
tstat_gamma_unfilt = (est_q_unfilt.b(3) - 10) / est_q_unfilt.se(3);
tstat_alpha_unfilt = (est_q_unfilt.b(1) - 0) / est_q_unfilt.se(1);
pval_gamma_unfilt = 1 - normcdf(tstat_gamma_unfilt);
pval_alpha_unfilt = 2 * (1 - normcdf(abs(tstat_alpha_unfilt)));

% Implied returns 
sdf_q_unfilt = beta_q * X_q_unfiltered(:,1) .^ (-est_q_unfilt.b(3));
ret_model_q_unfilt = est_q_unfilt.b(1) - mean((sdf_q_unfilt - est_q_unfilt.b(2)) .* r_excess_q, 1) ./ est_q_unfilt.b(2);
ret_model_pct_q_unfilt = ret_model_q_unfilt * 100;                     % to get teh percentage points

% Plot Unfiltered 
figure;
scatter(ret_model_pct_q_unfilt, ret_sample_pct_q, 30, [1 0 1], 'filled');  % magenta dots
hold on;
diag_min = min([ret_model_pct_q_unfilt, ret_sample_pct_q]);
diag_max = max([ret_model_pct_q_unfilt, ret_sample_pct_q]);
plot([diag_min diag_max], [diag_min diag_max], 'k-', 'LineWidth', 0.5);  % solid thin black line
xlabel('Model-Implied Excess Return (%)');
ylabel('Realized Mean Excess Return (%)');
title('Unfiltered Consumption: Model vs. Data');
xlim([diag_min diag_max]); ylim([diag_min diag_max]);
box on;
hold off;

% Interpretation of plot : Please refer pdf

% R square 
r2_q_unfilt = 1 - var(ret_sample_q - ret_model_q_unfilt) / var(ret_sample_q);

fprintf('\n==============================================================\n');
fprintf('        GMM Estimation Summary: Quarterly CBM Results         \n');
fprintf('==============================================================\n');

%% ------------------- Filtered Consumption -------------------
fprintf('\n>> Using Filtered Consumption Growth:\n');
disp(results_q_filt);
fprintf('  t-stat (gamma - 10):        %.4f\n', tstat_gamma_filt);
fprintf('  p-value (H0: gamma ≤ 10):   %.4f\n', pval_gamma_filt);
fprintf('  t-stat (alpha = 0):         %.4f\n', tstat_alpha_filt);
fprintf('  p-value (H0: alpha = 0):    %.4f\n', pval_alpha_filt);
fprintf('  Cross-sectional R²:         %.4f\n', r2_q_filt);
fprintf('  J-statistic:                %.4f\n', est_q_filt.J);

%% ------------------ Unfiltered Consumption ------------------
fprintf('\n>> Using Unfiltered Consumption Growth:\n');
disp(results_q_unfilt);
fprintf('  t-stat (gamma - 10):        %.4f\n', tstat_gamma_unfilt);
fprintf('  p-value (H0: gamma ≤ 10):   %.4f\n', pval_gamma_unfilt);
fprintf('  t-stat (alpha = 0):         %.4f\n', tstat_alpha_unfilt);
fprintf('  p-value (H0: alpha = 0):    %.4f\n', pval_alpha_unfilt);
fprintf('  Cross-sectional R²:         %.4f\n', r2_q_unfilt);
fprintf('  J-statistic:                %.4f\n', est_q_unfilt.J);

fprintf('\n==============================================================\n');
fprintf('         Interpretation and plots included in PDF file        \n');
fprintf('==============================================================\n');
