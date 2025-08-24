function [gT, uT] = Q_cbm_moments(params, ~, ~, Y, X, Z)
% -------------------------------------------------------------------------
% Q_cbm_moments
% Computes GMM moment vector and residual matrix for Quarterly CBM
%
% INPUTS:
%   params : 3x1 parameter vector [alpha; mu; gamma]
%   X      : T x (1+N) matrix = [cons_growth, excess_returns]
%   Y, Z   : not used (set to ~), included for toolbox compatibility
%
% OUTPUTS:
%   gT : (N+1)x1 vector of average moment conditions
%   uT : T x (N+1) matrix of residuals for each time period
% -------------------------------------------------------------------------

% Fixed model parameter
beta = 0.95;

% Extract parameters
alpha = params(1);    % intercept (alpha)
mu    = params(2);    % mean of the SDF
gamma = params(3);    % risk aversion coefficient

% Split data
cons_growth = X(:,1);             % consumption growth
asset_returns = X(:,2:end);       % excess returns for 25 portfolios

% Compute SDF based on power utility
sdf = beta * cons_growth .^ (-gamma);  % T x 1 vector

% Moment conditions:
%   1. Pricing errors: asset_returns - alpha + ((m_t - mu) * returns) / mu
pricing_errors = asset_returns - alpha + ((sdf - mu) .* asset_returns) / mu;

%   2. Mean-zero constraint on the SDF: E[m_t - mu] = 0
sdf_mean_constraint = sdf - mu;

% Stack all residuals
uT = [pricing_errors, sdf_mean_constraint];   % T x (25+1)

% Compute mean moment condition (to be minimized)
gT = mean(uT, 1)';    % (N+1)x1 column vector

end
