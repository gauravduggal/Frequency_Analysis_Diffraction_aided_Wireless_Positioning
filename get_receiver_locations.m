function receiver_locations = get_receiver_locations(data)
    % Initialize an empty container for receiver locations
    receiver_locations = zeros(3,length(data));

    % Loop through each receiver point in the data
    for i = 1:length(data)
        receiver_point = data(i).ReceiverPoint;
        
        % Check if the receiver has any paths
        paths = data(i).Paths;
        
        if ~isempty(paths)
            % Get the last interaction point from the first path
            last_interaction_point = paths(1).InteractionPoints(end, :);
            
            % Store the last interaction point as the receiver's location
            receiver_locations(:,i) = last_interaction_point.';
        else
            fprintf('Receiver point %d has no paths.\n', receiver_point);
        end
    end
end