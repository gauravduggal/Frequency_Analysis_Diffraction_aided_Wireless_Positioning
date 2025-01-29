% Function to calculate the first receiver path lengths
function [first_path] = get_first_receiver_path(rxdata, fap_threshold)
% receiver_path_lengths = containers.Map();  % Initialize an empty map
% receiver_path_lengths = zeros(1,length(data));
c = physconst('LightSpeed');
% interaction_str = cell(1,length(data));
% for i = 1:length(data)
receiver_point = rxdata.ReceiverPoint;

if isempty(rxdata.Paths)
    fprintf('Receiver %d has no paths.\n', receiver_point);
    return;
end

% Get the path with the lowest TOA (first path)
paths = rxdata.Paths;
data = [[paths.TimeOfArrival_sec_]',[paths.ReceivedPower_dBm_]'];

% Set the threshold for column 2
threshold = max(data(:,2))-fap_threshold;

% Find indices where column 2 is above the threshold
indices = find(data(:, 2) >= threshold);

% Filter rows using the indices
filtered_rows = data(indices, :);

% Check if there are any rows meeting the criteria
if ~isempty(filtered_rows)
    % Find the index of the minimum value in column 1 within filtered rows
    [~, minIndex] = min(filtered_rows(:, 1));

    % Get the original index in data
    originalIndex = indices(minIndex);

    
    % Display the result with the original index
    % disp('Row with smallest value in column 1 and column 2 above threshold:');
    % fprintf('Original Index: %d\n', originalIndex);
    % disp(result);
    % disp('No rows found with column 2 above the threshold.');
end

% [~, min_idx] = min([paths.TimeOfArrival_sec_]);
first_path = paths(originalIndex);

% Store the total path length for this receiver point
% receiver_path_lengths(i) = first_path.TimeOfArrival_sec_*c;
% interaction_str{i} = first_path.InteractionSummary;
% end
end