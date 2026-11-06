%% estimate_digital_concentration.m
% By Rena Fukuda (rf384@cornell.edu) for the manuscript 
% "A digital droplet microarray for measuring tolerance to antibiotics"

% Based on the framework in "Theoretical design and analysis of multivolume
% digital assays with wide dynamic range validated experimentally with 
% microfluidic digital PCR. "

% Purpose: Estimate four sample concentration (lambdas) from a digital 
% (positive/negative) SPOTs assay run at several volumes.
%
% Model
%   A volume v is negative with probability exp(-lambda*v).
%   lambda is found per sample by maximum likelihood (Newton's method with
%   backtracking), using the positive/negative counts at every volume.
%
% Input (Excel sheet)
%   - A row of well volumes (uL), one value per column.
%   - A digital matrix with the same columns: 1 = positive, 0 = negative,
%     NaN/empty = no well. Each sample occupies `rows_per_sample` rows.
%
% Output (workspace)
%   lambda  num_samples x 1   concentration, per mL (volumes are converted
%                             from uL to mL). -Inf = no positive wells,
%                             Inf = no negative wells.
%   sigma   num_samples x 1   approximate standard error of log(lambda),
%                             from the Fisher information. Inf if lambda
%                             is not estimable.
%
% Note: this code was cleaned by Claude Chat in November 2026.
clear;clc;

%% Settings
filename        = "26_02_20_Timecourse.xlsx";
sheet           = 13;
volume_range    = "A95:DJ95";   % well volume (uL) for each column
digital_range   = "A78:DJ93";   % 1 = positive, 0 = negative
rows_per_sample = 4;            % rows in the digital matrix per sample

% Newton solver settings
opts.max_iter = 100;            % max Newton iterations
opts.grad_tol = 1e-12;          % converged when |gradient| < grad_tol
opts.step_tol = 1e-10;          % converged when |Newton step| < step_tol

%% Load data
volume_row = readtable(filename, 'Sheet', sheet, 'Range', volume_range, ...
    'TreatAsEmpty', 'NaN');
volume_row = volume_row{1, :};

digital = table2array(readtable(filename, 'Sheet', sheet, ...
    'Range', digital_range));

assert(size(digital, 2) == numel(volume_row), ...
    "Volume row and digital matrix must have the same number of columns.");
assert(mod(size(digital, 1), rows_per_sample) == 0, ...
    "Digital matrix rows must be a multiple of rows_per_sample.");

%% Count positive and negative wells per sample and volume
vol_list    = unique(volume_row(~isnan(volume_row)));   % distinct volumes (uL)
num_vols    = numel(vol_list);
num_samples = size(digital, 1) / rows_per_sample;

pos_wells = NaN(num_samples, num_vols);
neg_wells = NaN(num_samples, num_vols);

for iv = 1:num_vols
    % Columns at this volume, regrouped as [rows_per_sample x samples x columns]
    cols   = volume_row == vol_list(iv);
    blocks = reshape(digital(:, cols), rows_per_sample, num_samples, []);

    pos_wells(:, iv) = reshape(sum(blocks == 1, [1 3]), [], 1);
    neg_wells(:, iv) = reshape(sum(blocks == 0, [1 3]), [], 1);
end

vols_ml = vol_list(:).' / 1000;   % volumes were entered in uL; convert to mL

%% Fit lambda for each sample
lambda = zeros(num_samples, 1);
sigma  = zeros(num_samples, 1);

for is = 1:num_samples
    [lambda(is), sigma(is)] = fit_lambda( ...
        pos_wells(is, :), neg_wells(is, :), vols_ml, opts);
end

%% Local functions
function [lambda, sigma] = fit_lambda(pos, neg, vols, opts)
% Maximum-likelihood lambda for one sample.
%   pos, neg : 1 x num_vols counts of positive / negative wells
%   vols     : 1 x num_vols well volumes (mL)
    tot = pos + neg;

    % Special cases: lambda is not estimable
    if sum(pos) == 0          % no positive wells
        lambda = -Inf;  sigma = Inf;
        return
    elseif sum(neg) == 0      % no negative wells
        lambda = Inf;   sigma = Inf;
        return
    end

    % Initial guess: lambda implied by each volume, clipped, then the max
    guess = log(1 - pos ./ tot) ./ -vols;
    guess(guess > 1e5)  = 1e5;
    guess(guess < 1e-5) = 1e-5;
    lambda = max(guess);

    % Newton iterations
    for iter = 1:opts.max_iter
        Ez = min(exp(-lambda * vols), 1 - 1e-10);   % clip for safety

        % Gradient and Hessian of the log-likelihood
        g = sum(-vols .* (tot - pos ./ (1 - Ez)));
        if abs(g) < opts.grad_tol, break, end

        h    = -sum(pos .* vols.^2 .* Ez ./ (1 - Ez).^2);
        step = -g / h;
        if abs(step) < opts.step_tol, break, end

        % Keep lambda non-negative
        while lambda + step < 0
            step = step / 2;
        end

        % Backtracking: shrink the step while the likelihood gets much worse
        f0 = loglik(lambda, pos, neg, vols);
        while loglik(lambda + step, pos, neg, vols) < f0 - 0.1*abs(f0) - 1e-5
            step = step / 2;
        end

        lambda = lambda + step;
    end

    % Uncertainty: Fisher information for log(lambda) (includes lambda^2 factor)
    Ez   = exp(-lambda * vols);
    info = sum(tot .* (lambda * vols).^2 .* Ez ./ (1 - Ez));

    if info <= 0
        sigma = Inf;
    else
        sigma = 1 / sqrt(info);
    end
end

function f = loglik(lambda, pos, neg, vols)
% Log-likelihood of lambda (up to a constant).
    Ez = min(exp(-lambda * vols), 1 - 1e-14);
    f  = sum(pos .* log(1 - Ez) + neg .* (-lambda * vols));
end