% Parameters
n_firms = 1000; % Number of firms
n_periods = 100; % Total periods
rho = 0.7; % AR(1) coefficient for productivity
sigma_eps = 0.3; % Std deviation of productivity shock
alphak = 0.2; alphal = 0.3; alpham = 0.5; % Cobb-Douglas parameters
rng(42); % For reproducibility

% Simulate capital (K) - lognormal distribution
K = lognrnd(0, 1, [n_firms, n_periods]);

% Simulate initial productivity (A_0)
A_0 = normrnd(0, 1, [n_firms, 1]);


% Generate productivity (A) using AR(1) process
A = zeros(n_firms, n_periods);
A(:, 1) = A_0;
for t = 2:n_periods
    A(:, t) = rho * A(:, t-1) + normrnd(0, sigma_eps, [n_firms, 1]);
end
A=exp(A);

L=(A.*alpham.^(1-alphal)*alphal.^alphal.*K.^alphak).^(1/(1-alpham-alphal));
M=(A.*alphal.^(1-alpham)*alpham.^alpham.*K.^alphak).^(1/(1-alpham-alphal));
Y = A .* (K .^ alphak) .* (L .^ alphal) .* (M .^ alpham);

% Extract the last 10 periods for analysis
periods_to_keep = (n_periods - 9):n_periods;
K_final = K(:, periods_to_keep);
A_final = A(:, periods_to_keep);
L_final = L(:, periods_to_keep);
M_final = M(:, periods_to_keep);
Y_final = Y(:, periods_to_keep);

% Assuming you have matrices for Capital, Productivity, Labor, Material, and Output
% from the previous simulation, each being n_firms x n_periods in size

% Log-transform the variables
LogCapital = log(K_final);
LogProductivity = log(A_final);
LogLabor = log(L_final);
LogMaterial = log(M_final);
LogOutput = log(Y_final);
% Flatten the matrices for correlation calculation since MATLAB corr function
% operates on vectors or matrices where each column represents a variable
LogCapitalFlat = LogCapital(:);
LogProductivityFlat = LogProductivity(:);
LogLaborFlat = LogLabor(:);
LogMaterialFlat = LogMaterial(:);
LogOutputFlat = LogOutput(:);

% Combine the log-transformed variables into a single matrix
dataMatrix = [LogCapitalFlat, LogProductivityFlat, LogLaborFlat, LogMaterialFlat, LogOutputFlat];

% Calculate the correlation matrix
correlationMatrix = corr(dataMatrix, 'Rows', 'complete');  % 'complete' to use rows with no NaNs

% Display the correlation matrix
correlationMatrix


% Assuming Output and Productivity matrices are n_firms x n_periods in size

% Step 1: Calculate Output Market Share (si) for each firm in each period
total_output_per_period = sum(Y_final, 1); % Sum across rows for each period
si = Y_final ./ total_output_per_period; % Broadcasting division to calculate market shares

% Step 2: For simplicity, let's use output share as weights (ωi = si)
% If using productivity or another metric, replace si with that metric here
wi = log(A_final); % This example uses si directly for simplicity

% Step 3: Compute Weighted Sum (ωisi) and Sum of Weights (ωi)
weighted_sum = sum(wi .* si, 'all'); % Multiplying element-wise and summing over all elements
sum_of_weights = sum(wi, 'all'); % Summing all weights

% Number of observations (N) is the total number of firm-period observations
N = numel(Y_final); % Total number of elements in Output matrix

% Step 4: Compute Olley-Pakes Covariance
%OP_covariance = weighted_sum - (sum_of_weights / N);


OP_covariance = (weighted_sum - sum(wi*1/N, 'all'))/10;
weighted_productivity = sum(sum(wi .* si)) / sum(sum(si)); % Omega
unweighted_productivity = mean(wi, 'all'); % omega_bar


% Display the Olley-Pakes Covariance
disp('Olley-Pakes Covariance:');
disp(OP_covariance);
disp('weighted_productivity:');
disp(weighted_productivity);
disp('unweighted_productivity:');
disp(unweighted_productivity);

