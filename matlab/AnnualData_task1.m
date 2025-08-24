clear          
clc            

seed = 1;      % fixed seed for reprodusibility
rng(seed);     % random number generator seed

%% Preprocessing
% Read both datasets (annual & quarterly)
AnnualTbl = readtable('AnnualData_1949_2018.csv');
QuarterlyTbl = readtable('QuarterlyData_1949Q1_2024Q1.csv');

% Extract year column and convert rest to numerical array
years_vec = AnnualTbl{:,1};
A_raw = table2array(AnnualTbl(:,2:end));    % size: 70x3

% Assign variables for clarity
cons_filtered = A_raw(:,1);                 % filtered consumption growth --> column 1
cons_unfiltered = A_raw(:,2);               % unfiltered consumptiongrowth --> column2 
ret_excess = A_raw(:,3);                    % excess return on market portfolio --> column 3

%% =========================================================================
%% SECTION 1.1 : Annual Data - Classical Consumption-Based Model (CBM)
%% ========================================================================

%% 1.1.1 --> A function that returns gT and uT is - A_cbm_momemts(Annual CBM moment)

%% 1.1.2 --> Outline the features of the Consumption-Based Model --> please refer PDF 

%% 1.1.3 Estimation of gamma using gmm toolbox 
addpath('gmm');
addpath('minz');

gmmopt = struct();
gmmopt.infoz.momt = 'A_cbm_moments';    % name of moment function
gmmopt.W0 = 'I';                        % identity matrix weighting
gmmopt.gmmit = 1;                       % one-step GMM
gmmopt.hess = 'gn';                    % Gauss-Newton
gmmopt.S = 'NW';                        % Newey-West spectral density
gmmopt.plot = 0;                        
gmmopt.prt = 1;                         

% Inputs
init = 3;                               % initial guess for gama
X = [cons_filtered, ret_excess];        % X = [c_growth, r_excess]
Y = zeros(70,1);                        % moment targets = 0
Z = ones(70,1);                         % constant instrument

% Run GMM
[est_filtered, res_filtered] = gmm(init, gmmopt, Y, X, Z);


% Format result table
tab_filtered = table(est_filtered.b, est_filtered.se, ...
    est_filtered.b - 1.96 * est_filtered.se, ...
    est_filtered.b + 1.96 * est_filtered.se, ...
    'VariableNames', {'Estimate', 'StdError', 'CI_Lower', 'CI_Upper'});

% disp(tab_filtered);

%% 1.1.4 Hypothesis Test: H0: gamma ≤ 10 (one-sided t-test)
t_stat_filtered = (est_filtered.b - 10) / est_filtered.se;                 % computes t-stat
p_value_filtered = 1 - normcdf(t_stat_filtered);                           % related p-value for one sided t-test

% disp(['t-statistic = ', num2str(t_stat_filtered)]);
% disp(['p-value (H0: gamma ≤ 10) = ', num2str(p_value_filtered)]);

%% 1.1.5 Compute Implied Risk-Free Rate from CBM
beta_bar = 0.95;
gamma_hat = est_filtered.b;
Rf_filtered = -log(beta_bar) + gamma_hat * mean(log(cons_filtered)) ...
              - 0.5 * gamma_hat^2 * var(log(cons_filtered));

% disp(['Implied Risk-Free Rate = ', num2str(Rf_filtered)]);

%% 1.1.6 Time Series of SDF and Visualization
SDF_filtered = beta_bar * cons_filtered.^(-gamma_hat);      % compute SDF values
sdf_dates = datetime(years_vec, 1, 1);                         

% Plot SDF over time with US recessions 
figure('Color', [0.95 0.95 0.95]);  % light grey figure background
plot(sdf_dates, SDF_filtered, 'Color', [1 0 1], 'LineWidth', 1.2);  
hold on;
ylim([0 11]);
recessionplot;
xlabel('Year', 'FontWeight', 'bold');
ylabel('SDF (m_t)', 'FontWeight', 'bold');
legend('SDF path', 'US recessions', 'Location', 'best');
title('SDF Path using Filtered Consumption', 'FontWeight', 'bold');
set(gca, 'Box', 'on', 'Color', [1 1 1]);  % white axes background
grid on;
hold off;


% Correlation with market return
corr_SDF_filtered = corr(SDF_filtered, ret_excess);

% disp(['Corr(SDF, Excess Return) = ', num2str(corr_SDF_filtered)]);

%% =========================================================================
%% 1.2 A New Measure for Consumption Growth 
%% =========================================================================

% Redefine X using unfiltered consumption measure
X_unfiltered = [cons_unfiltered, ret_excess];

% Run GMM estimation again with consistent naming
[est_unfiltered, res_unfiltered] = gmm(init, gmmopt, Y, X_unfiltered, Z);


% Result Table
tab_unfiltered = table(est_unfiltered.b, est_unfiltered.se, ...
    est_unfiltered.b - 1.96 * est_unfiltered.se, ...
    est_unfiltered.b + 1.96 * est_unfiltered.se, ...
    'VariableNames', {'Estimate', 'StdError', 'CI_Lower', 'CI_Upper'});

% disp(tab_unfiltered);

% Hypothesis Test: H0: gamma ≤ 10 (Unfiltered Consumption)
t_stat_unfiltered = (est_unfiltered.b - 10) / est_unfiltered.se;
p_value_unfiltered = 1 - normcdf(t_stat_unfiltered);

% disp(['t-statistic (unfiltered) = ', num2str(t_stat_unfiltered)]);
% disp(['p-value (H0: gamma ≤ 10, unfiltered) = ', num2str(p_value_unfiltered)]);

% Compute implied Risk-Free Rate 
gamma_hat_unfilt = est_unfiltered.b;
Rf_unfiltered = -log(beta_bar) + gamma_hat_unfilt * mean(log(cons_unfiltered)) ...
              - 0.5 * gamma_hat_unfilt^2 * var(log(cons_unfiltered));

% disp(['Implied Risk-Free Rate (unfiltered) = ', num2str(Rf_unfiltered)]);

% Estimate SDF over time 
SDF_unfiltered = beta_bar * cons_unfiltered.^(-gamma_hat_unfilt);

% Plot SDF Unfiltered 
figure('Color', [0.95 0.95 0.95]);  % light grey figure background
plot(sdf_dates, SDF_unfiltered, 'Color', [1 0 1], 'LineWidth', 1.2);  % magenta line
hold on;
ylim([0 11]);
recessionplot;
xlabel('Year', 'FontWeight', 'bold');
ylabel('SDF (m_t)', 'FontWeight', 'bold');
legend('SDF path', 'US recessions', 'Location', 'best');
title('SDF Path using Unfiltered Consumption', 'FontWeight', 'bold');
set(gca, 'Box', 'on', 'Color', [1 1 1]);  % white axes background
grid on;
hold off;

% Correlation with market return (Unfiltered)
corr_SDF_unfiltered = corr(SDF_unfiltered, ret_excess);

% disp(['Corr(SDF, Excess Return) — Unfiltered: ', num2str(corr_SDF_unfiltered)]);


%% =========================================================================
%% Final Display Block: Summary of Results (for Appendix)
%% =========================================================================

fprintf('\n======================== ANNUAL CBM ESTIMATION SUMMARY ========================\n');

% ---------- Filtered Consumption ----------
fprintf('\n[1] Filtered Consumption Growth:\n');
fprintf('Estimated gamma (risk aversion):       %.4f\n', est_filtered.b);
fprintf('Standard Error:                        %.4f\n', est_filtered.se);
fprintf('95%% Confidence Interval:               [%.4f, %.4f]\n', ...
    est_filtered.b - 1.96 * est_filtered.se, est_filtered.b + 1.96 * est_filtered.se);
fprintf('t-statistic (gamma - 10):              %.4f\n', t_stat_filtered);
fprintf('p-value (H0: gamma ≤ 10):              %.4f\n', p_value_filtered);
fprintf('Implied Risk-Free Rate:                %.4f\n', Rf_filtered);
fprintf('Corr(SDF, Excess Return):              %.4f\n', corr_SDF_filtered);

% ---------- Unfiltered Consumption ----------
fprintf('\n[2] Unfiltered Consumption Growth:\n');
fprintf('Estimated gamma (risk aversion):       %.4f\n', est_unfiltered.b);
fprintf('Standard Error:                        %.4f\n', est_unfiltered.se);
fprintf('95%% Confidence Interval:               [%.4f, %.4f]\n', ...
    est_unfiltered.b - 1.96 * est_unfiltered.se, est_unfiltered.b + 1.96 * est_unfiltered.se);
fprintf('t-statistic (gamma - 10):              %.4f\n', t_stat_unfiltered);
fprintf('p-value (H0: gamma ≤ 10):              %.4f\n', p_value_unfiltered);
fprintf('Implied Risk-Free Rate:                %.4f\n', Rf_unfiltered);
fprintf('Corr(SDF, Excess Return):              %.4f\n', corr_SDF_unfiltered);

% ---------- Conclusion ----------
fprintf('\n[3] Summary Comment:\n');
if corr_SDF_unfiltered > corr_SDF_filtered
    disp('→ Unfiltered consumption co-varies better with excess returns.');
else
    disp('→ Filtered consumption shows stronger SDF-return comovement.');
end

fprintf('==============================================================================\n');
