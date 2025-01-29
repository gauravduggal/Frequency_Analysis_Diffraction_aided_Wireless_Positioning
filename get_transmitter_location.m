function transmitter_location = get_transmitter_location(data)
    % Initialize transmitter location
    transmitter_location = [];

    % Check if the first receiver point exists
    if ~isempty(data)
        % Get the first receiver point
        first_receiver = data(1);
        
        % Check if it has paths
        paths = first_receiver.Paths;
        if ~isempty(paths)
            % Get the first interaction point from the first path (assumed to be the transmitter)
            transmitter_location = paths(1).InteractionPoints(1, :).';
        else
            fprintf('No paths found for the first receiver point.\n');
        end
    else
        fprintf('No receiver points found in data.\n');
    end
end
