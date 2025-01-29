% Function to calculate the first receiver path lengths
function [mpc_group_FAP] = get_mpc_group_FAP(rxdata, mpc_group)
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
Np = length(paths);

flag_found = false;
min_tof = 0;
min_tof_idx = 0;
for pidx = 2:Np
    for tyidx = 1:length(mpc_group)
    if strcmpi(paths(pidx).InteractionSummary, mpc_group{tyidx})
       if flag_found==false
           min_tof = paths(pidx).TimeOfArrival_sec_;
           min_tof_idx = pidx;
           flag_found=true;
       else
           if min_tof > paths(pidx).TimeOfArrival_sec_
               min_tof = paths(pidx).TimeOfArrival_sec_;
               min_tof_idx = pidx;
           end
       end
    end
    end
end
if min_tof_idx~=0
    mpc_group_FAP = paths(min_tof_idx);
else
    mpc_group_FAP = paths(1);
    fprintf('Receiver %d has no diffraction paths.\n', receiver_point);
    return
end

return;




% Store the total path length for this receiver point
% receiver_path_lengths(i) = first_path.TimeOfArrival_sec_*c;
% interaction_str{i} = first_path.InteractionSummary;
% end
end