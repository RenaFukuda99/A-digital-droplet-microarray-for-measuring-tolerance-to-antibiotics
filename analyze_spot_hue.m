%% analyze_spot_hue.m
% By Rena Fukuda (rf384@cornell.edu) for the manuscript 
% "A digital droplet microarray for measuring tolerance to antibiotics"
% Note: code has been cleaned by Claude Chat in November 2026 

% Ojective: Extract the mean hue of every spot on a plate image.
%
% Workflow
%   1. Click 18 calibration points on the image (9 per spot set).
%   2. Spot centers are interpolated between the clicked points.
%   3. The mean HSV hue inside a circle around each center is computed.
%   4. The spot sets are interleaved back into the plate layout and
%      written to a spreadsheet.
%
% Click order (set 1 first, then set 2; press Enter when finished)
%   For each set, click the Left, Middle, and Right spot of the
%   top row -> middle row -> bottom row (9 clicks per set, 18 total).
%
% Requires: Image Processing Toolbox (getpts, viscircles)

clc; clear; close all;

%% Settings
image_file      = "myfile.jpeg";
output_file     = "output.xlsx";

num_rows        = 4;     % rows per spot set (must be >= 3)
num_cols_set1   = 29;    % columns in spot set 1
num_cols_set2   = 28;    % columns in spot set 2

% Which set supplies the first column of the output (1 or 2).
% The first set must have the same number of columns as the other set,
% or exactly one more.
first_in_output = 1;

sample_radius   = 10;    % px, radius of the area averaged for each spot
display_radius  = 5;     % px, radius of the circles drawn as a visual check

%% Load image and collect calibration points
im_rgb = imread(image_file);
hue    = rgb2hsv(im_rgb);
hue    = hue(:,:,1);

imshow(im_rgb, 'InitialMagnification', 'fit');
title("Click 9 points for set 1, then 9 for set 2, then press Enter");

[x, y] = getpts;
assert(numel(x) == 18, "Expected 18 clicked points, got %d.", numel(x));

%% Interpolate spot centers
x_set1 = interpolate_grid(x(1:9),   num_rows, num_cols_set1);
y_set1 = interpolate_grid(y(1:9),   num_rows, num_cols_set1);
x_set2 = interpolate_grid(x(10:18), num_rows, num_cols_set2);
y_set2 = interpolate_grid(y(10:18), num_rows, num_cols_set2);

% All centers as [x y]; set 1 first, each set in column-major order
centers   = [x_set1(:), y_set1(:); x_set2(:), y_set2(:)];
num_spots = size(centers, 1);

% Visual check of spot positions
hold on
viscircles(centers, display_radius);

%% Build a label image (pixel value = spot index, 0 = background)
[xx, yy]   = meshgrid(1:size(hue, 2), 1:size(hue, 1));
spot_label = zeros(size(hue), 'uint16');

for ii = 1:num_spots
    in_spot = hypot(xx - centers(ii,1), yy - centers(ii,2)) <= sample_radius;
    spot_label(in_spot) = ii;
end

% Show the regions that will be sampled
figure;
imshow(im_rgb .* cast(spot_label > 0, 'like', im_rgb), 'InitialMagnification', 'fit');
title("Sampled regions");

%% Mean hue per spot
mean_hue = nan(1, num_spots);
for ii = 1:num_spots
    mean_hue(ii) = mean(hue(spot_label == ii));
end

%% Reassemble into the plate layout
n_set1   = num_rows * num_cols_set1;
hue_set1 = reshape(mean_hue(1:n_set1),     num_rows, num_cols_set1);
hue_set2 = reshape(mean_hue(n_set1+1:end), num_rows, num_cols_set2);

if first_in_output == 1
    hue_plate = interleave_columns(hue_set1, hue_set2);
else
    hue_plate = interleave_columns(hue_set2, hue_set1);
end

%% Save results
writematrix(hue_plate, output_file);
fprintf("Saved %d x %d hue matrix to %s\n", ...
    size(hue_plate, 1), size(hue_plate, 2), output_file);

%% Local functions
function grid = interpolate_grid(p, num_rows, num_cols)
% Interpolate a num_rows x num_cols grid of coordinates (x or y) from 9
% clicked values, ordered: top L/M/R, middle L/M/R, bottom L/M/R.
    assert(num_rows >= 3 && num_cols >= 3, ...
        "Grid needs at least 3 rows and 3 columns.");

    half_col = floor(num_cols / 2);
    half_row = floor(num_rows / 2);
    grid     = zeros(num_rows, num_cols);

    % Interpolate along the top, middle, and bottom rows
    row_idx = [1, half_row, num_rows];
    for k = 1:3
        left  = p(3*k - 2);
        mid   = p(3*k - 1);
        right = p(3*k);
        r = row_idx(k);
        grid(r, 1:half_col)   = linspace(left, mid,   half_col);
        grid(r, half_col:end) = linspace(mid,  right, num_cols - half_col + 1);
    end

    % Interpolate down each column (top -> middle -> bottom)
    for c = 1:num_cols
        grid(1:half_row, c) = linspace(grid(1, c), grid(half_row, c), half_row);
        grid(half_row:num_rows, c) = linspace(grid(half_row, c), ...
            grid(num_rows, c), num_rows - half_row + 1);
    end
end

function out = interleave_columns(A, B)
% Interleave columns as A1 B1 A2 B2 ...
% A must have the same number of columns as B, or one more.
    nB = size(B, 2);
    assert(size(A, 2) == nB || size(A, 2) == nB + 1, ...
        "First set must have equal or one more column than the second.");

    paired = reshape([A(:, 1:nB); B], size(A, 1), []);
    out    = [paired, A(:, nB+1:end)];   % trailing column of A, if any
end