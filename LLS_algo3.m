function [np_est,res] = LLS_algo3(r,a,np)
% Linear Least Squares 3D Positioning with All Measurement Differences

anchors=a';
% Number of anchors
Na = size(anchors, 1);

% Initialize matrices for pairwise differences
num_pairs = nchoosek(Na, 2); % Number of unique anchor pairs
A = zeros(num_pairs, 3); % Design matrix
b = zeros(num_pairs, 1); % Difference vector

% Generate pairwise equations
row = 1; % Row index for A and b
for i = 1:Na-1
    for j = i+1:Na
        % Compute the row for A and b
        A(row, :) = 2 * (anchors(i, :) - anchors(j, :));
        b(row) = r(j)^2 - r(i)^2 ...
                 - sum(anchors(j, :).^2) + sum(anchors(i, :).^2);
        row = row + 1;
    end
end

% Solve for the target position using linear least squares
np_est = (A \ b);

res = sqrt(sum((np_est-np).^2));
end