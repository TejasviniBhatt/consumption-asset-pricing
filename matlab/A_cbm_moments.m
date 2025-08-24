function [gT, uT] = A_cbm_moments(b, ~, ~, Y, X, Z)
% -------------------------------------------------------------------------
% A_cbm_moments
% Computes the moment condition and residuals for CBM using Annual Data
% Inputs:
%   b : scalar parameter (γ̃)
%   X : [T x 2] matrix with consumption growth and excess returns
% Outputs:
%   gT : mean moment condition (scalar)
%   uT : residuals for each time point (T x 1)
% -------------------------------------------------------------------------

% Model parameters
beta = 0.95;            % fixed subjective discount factor
gamma = b;              % input parameter: risk aversion coefficient

% Parse input data
c_growth = X(:,1);      % column 1: consumption growth
r_excess = X(:,2);      % column 2: market excess returns

% Compute time-series of moment conditions (residuals)
uT = (beta * c_growth.^(-gamma)) .* r_excess;

% Compute average moment (GMM criterion)
gT = mean(uT);

end
