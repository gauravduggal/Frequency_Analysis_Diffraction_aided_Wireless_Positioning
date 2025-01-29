clc
clear
close all
root_path = 'Data\';
output_file_path = 'Results\';
%%
% Refer to https://arxiv.org/pdf/2409.02832 for Diffraction based
% Wireless positioning technique


%%
if ~exist(output_file_path, 'dir')
    mkdir(output_file_path); % Create the directory
    fprintf('Directory "%s" created.\n', output_file_path);
else
    fprintf('Directory "%s" already exists.\n', output_file_path);
end

%% Loading data
frequencies = {'0.8 GHz', '2.4 GHz', '5.8 GHz','8 GHz', '10 GHz', '15 GHz', '28 GHz','37 GHz','48 GHz'};
colors = {
    [0.0000, 0.0000, 1.0000];  % Bright Blue
    [1.0000, 0.0000, 0.0000];  % Bright Red
    [0.0000, 1.0000, 0.0000];  % Bright Green
    [0.0000, 0.0000, 0.0000];  % Black
    [1.0000, 0.0000, 1.0000];  % Magenta
    [0.0000, 1.0000, 1.0000];  % Cyan
    [0.5000, 0.0000, 1.0000];  % Purple
    [1.0000, 0.5000, 0.0000];  % Orange
    [0.5000, 0.5000, 0.5000];  % Gray
    };
fap_threshold_vec = [0,10,20];
Nthresh = length(fap_threshold_vec);
mpc3 = {'Tx-D-Rx','Tx-D-T-Rx','Tx-D-T-T-Rx','Tx-D-T-T-T-Rx','Tx-D-T-T-T-T-Rx'};


base_path_vec = {
    '\Ray Tracing Simulation Code_800MHz_working\0p8GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_2p4GHz_working\2p4GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_5p8GHz_working\5p8GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_2R_2T_2D_APG_2m',
    '\Ray Tracing Simulation Code_8GHz_working\8GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_10GHz_working\10GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_15GHz_working\15GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_28GHz_working\28GHz_4TX_30dBm_Combined_Data_4Anchors_diffused_mpc_2R_2T_2D',
    '\Ray Tracing Simulation Code_37GHz_working\37GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_48GHz_working\48GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m'
    };
Nf = length(base_path_vec);
%number of anchors
Na = 4;
base_path = base_path_vec{1};
load(sprintf('%s%s%s',root_path,base_path,'\Transmitter1_Combined.mat'),'data');
tx_data(:,1) = data;
Nr = size(tx_data,1);

%% Simulation Parameters
w = 2;
%Span of Building

x1 = -15;
x2 = 15;
%Bandwidth
B = 400e6;
rx_processing_gain_vec=[-10,-10,-10,0,0,0,20,20,20];

k = physconst('Boltzmann');
T = 300;
N0 = k*T*B;
c = physconst('LightSpeed');

pebliml = 0;
peblimh = 15;
dim=1:3;
edges = pebliml:0.005:peblimh;
Nbins = length(edges)-1;


%%
data1 = zeros(Nbins,Nf,Nthresh);
data2 = zeros(Nbins,Nf,Nthresh);
data3 = zeros(Nbins,Nf,Nthresh);

for base_path_idx = 1:Nf

    rx_processing_gain = rx_processing_gain_vec(base_path_idx);

    base_path = base_path_vec{base_path_idx};
    load(sprintf('%s%s%s',root_path,base_path,'\Transmitter1_Combined.mat'),'data');
    tx_data(:,1) = data;
    load(sprintf('%s%s%s',root_path,base_path,'\Transmitter2_Combined.mat'), 'data');
    tx_data(:,2) = data;
    load(sprintf('%s%s%s',root_path,base_path,'\Transmitter3_Combined.mat'), 'data');
    tx_data(:,3) = data;
    if Na==4
        load(sprintf('%s%s%s',root_path,base_path,'\Transmitter4_Combined.mat'), 'data');
        tx_data(:,4) = data;
    end
    rx_loc = get_receiver_locations(tx_data(:,1));

    clear data
    % Anchor positions
    tx_loc(:,1) = get_transmitter_location(tx_data(:,1));
    tx_loc(:,2) = get_transmitter_location(tx_data(:,2));
    tx_loc(:,3) = get_transmitter_location(tx_data(:,3));
    if Na==4
        tx_loc(:,4) = get_transmitter_location(tx_data(:,4));
    end

    for fap_threshold_idx = 1:Nthresh
        fap_threshold = fap_threshold_vec(fap_threshold_idx);


        %For rx locations skipped due to low snr out a very high error
        peb_diff = (peblimh+1)*ones(1,Nr);
        rmse_ippa = (peblimh+1)* ones(1,Nr);
        rmse_nls = (peblimh+1)* ones(1,Nr);
        rmse_lls = (peblimh+1)* ones(1,Nr);
        rmse_nls_delta = (peblimh+1)* ones(1,Nr);




        %anchor positions
        a = zeros(size(tx_loc));
        %anchor coordinates transformed from WI to TWC paper model
        %y coordinate (LHS)
        a(2,:) = tx_loc(1,:);
        %x coordinate (LHS)
        a(1,:) = tx_loc(2,:);
        %z coordinate (LHS)
        a(3,:) = tx_loc(3,:);


        parfor tidx = 1:Nr
            flag_rx_skip = false;
            %Convert Wireless Insight to manuscript coordinate system
            xn = rx_loc(2,tidx);
            yn = rx_loc(1,tidx);
            zn = rx_loc(3,tidx);
            np = [xn;yn;zn];
            %we consider the lower diffraction MPC in our model
            zb = zn-w/2;

            FIM_diff = zeros(3,3);
            r = zeros(Na,1);
            weights = zeros(Na,Na);

            for aidx=1:Na
                %pick the first arriving path
                if tx_data(tidx,aidx).NumberOfPaths==0
                    continue
                end
                FAP = get_first_receiver_path(tx_data(tidx,aidx),fap_threshold);
                mpc3_FAP = get_mpc_group_FAP(tx_data(tidx,aidx),mpc3);


                rxPowerdBm = FAP.ReceivedPower_dBm_ + rx_processing_gain;
                mpc3_FAPrxPowerdbM = mpc3_FAP.ReceivedPower_dBm_ + rx_processing_gain;

                rxPower = 10^((rxPowerdBm)/10)*1e-3;
                mpc3_FAP_rxPower = 10^((mpc3_FAPrxPowerdbM)/10)*1e-3;

                snr_linear = rxPower/N0;
                mpc3_FAP_snr_linear = mpc3_FAP_rxPower/N0;
                weights(aidx,aidx) = snr_linear;
                %if SNR is very low, FIM can become non invertible due to numerical
                %errors so we assume RMSE is 100 at that location
                % if 10*log10(snr_linear) < -20
                %     flag_rx_skip = true;
                %     break
                % end

                xa = tx_loc(2,aidx);
                ya = tx_loc(1,aidx);
                za = tx_loc(3,aidx);

                [Qe,p,dp_dxn,dp_dyn,dp_dzn] = get_diffraction_coord_fermat(xa,ya,za,x1,x2,xn,yn,zn,zb,0.01);
                J_xx = dp_dxn^2;
                J_yy = dp_dyn^2;
                J_zz = dp_dzn^2;

                J_xy = dp_dxn*dp_dyn;
                J_yx = J_xy;

                J_yz = dp_dyn*dp_dzn;
                J_zy = J_yz;

                J_xz = dp_dxn*dp_dzn;
                J_zx = J_xz;

                J = [J_xx,J_yx,J_zx;
                    J_xy,J_yy,J_zy;
                    J_xz,J_yz,J_zz;];
                k = 8*pi^2*mpc3_FAP_snr_linear*B^2/c^2;
                FIM_diff = FIM_diff + k*J;
                % trace(inv(FIM_diff))

                % p = sqrt((xa-xn)^2+(ya-yn)^2+(za-zn)^2);
                % J_xx = (xa-xn)^2/p^2;
                % J_yy = (ya-yn)^2/p^2;
                % J_zz = (za-zn)^2/p^2;
                % J_xy = (xa-xn)*(ya-yn)/p^2;
                % J_yz = (ya-yn)*(za-zn)/p^2;
                % J_xz = (xa-xn)*(za-zn)/p^2;
                % J_zx = J_xz;
                % J_yx = J_xy;
                % J = [J_xx,J_yx,J_zx;
                %     J_xy,J_yy,J_zy;
                %     J_xz,J_yz,J_zz;];
                % k = 8*pi^2*snr_linear*B^2/c^2;
                % FIM_los = FIM_los + k*J;
                % r(aidx) = p + randn(1)/sqrt(k/2);
                r(aidx) = c*FAP.TimeOfArrival_sec_ + randn(1)/sqrt(k/2);


            end
            weights = weights/trace(weights);
            % if flag_rx_skip ==true
            %     continue
            % end
            % initial_estimate = [0;20;5];
            % initial_estimate = 100*rand(3,1);
            [np_est_lls,~] = LLS_algo3(r,a,np);

            rmse_lls_temp = sqrt(sum((np_est_lls(dim)-np(dim)).^2));

            % initial_estimate = np_est_lls;
            initial_estimate = [0;10;15];
            % initial_estimate = [30;20;5].*(rand(3,1)-[0.5;0;0])+[0;0;10];
            % initial_estimate(1) = 0;

            [np_est_nls1, res1] = nls_3D_estimator(r,a,initial_estimate, 10, w, x1,x2,weights,np);
            for sidx=1:100
                % initial_estimate = [30;20;5].*(rand(3,1)-[0.5;0;0])+[0;0;10];
                initial_estimate = np_est_nls1 + 10*(rand(3,1)-[0.5;0.5;0.5]);
                [np_est_nls, res] = nls_3D_estimator(r,a,initial_estimate, 10, w, x1,x2,weights,np);
                if res<res1
                    res1 = res;
                else
                    np_est_nls = np_est_nls1;
                    res = res1;
                end
            end
            rmse_nls_temp = sqrt(sum((np_est_nls(dim)-np(dim)).^2));



            A = FIM_diff; % Replace with your matrix
            conditionNumber = cond(A);

            if conditionNumber > 1e10 % or another threshold based on your requirement
                disp('Matrix is close to being singular (ill-conditioned)');
            end
            FIM_diff_inv = inv(FIM_diff);
            temp = diag(FIM_diff_inv);


            peb_diff(1,tidx) = abs(sqrt(sum(temp(dim))));%sqrt(FIM_diff_inv(1,1)+FIM_diff_inv(2,2)+FIM_diff_inv(3,3));
            rmse_lls(1,tidx) = rmse_lls_temp;
            rmse_nls(1,tidx) = rmse_nls_temp;

        end

        val1 = peb_diff(1,:);
        h1 = histogram(val1(:),edges,'Normalization','cdf','DisplayStyle','stairs','LineStyle','-.','LineWidth',2);
        data1(:,base_path_idx,fap_threshold_idx) = h1.Values;

        val2 = rmse_nls(1,:);
        h2 = histogram(val2(:),edges,'Normalization','cdf','DisplayStyle','stairs','LineStyle','-.','LineWidth',2);
        data2(:,base_path_idx,fap_threshold_idx) = h2.Values;


        val3 = rmse_lls(1,:);
        h3 = histogram(val3(:),edges,'Normalization','cdf','DisplayStyle','stairs','LineStyle','-.','LineWidth',2);
        data3(:,base_path_idx,fap_threshold_idx) = h3.Values;


        close all

    str = sprintf("%s:threshold:%d",base_path, fap_threshold);
    disp(str)
    end
end

%% FR1
% Define marker vector for unique plot markers
markervec = {'o', 'diamond', 'x', '^', 'pentagram'};


for fap_threshold_idx=1:Nthresh
    fap_threshold = fap_threshold_vec(fap_threshold_idx);
    % Create figure with specific dimensions for publication
    figure('Units', 'centimeters' , 'Position', [1 1 16.8 12], 'PaperPositionMode', 'auto');


    markerindices = 1:200:Nbins;


    for fileidx=1:3
        % Plot data with specified properties
        legend_str = sprintf("CRLB:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data1(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{1}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx});

        hold on;
        legend_str = sprintf("D-NLS:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data2(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{2}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx}); % Red
        legend_str = sprintf("LLS:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data3(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{5}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx}); % Yellow



        % Customize axes
        xlim([0, edges(end)]);
        ylim([0, 1]);
        grid on;
        xlabel('RMSE (m)', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('CDF', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off');

        ax = gca; % Get the current axes
        ax.XMinorGrid = 'on'; % Enable minor grid for x-axis
        ax.YMinorGrid = 'on'; % Enable minor grid for y-axis

        % Set the spacing for major ticks
        ax.XTick = 0:2:peblimh;
        ax.YTick = 0:0.1:1;

        % Aspect ratio for the plot
        % daspect([5 0.8 1]);
        daspect([10 0.8 1]);

    end
    % Customize the legend for clarity
    % legend('Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 14);
    legend('Interpreter', 'latex', 'Location', 'northoutside', 'Orientation', 'horizontal','NumColumns', 3,'FontSize', 10,'Box', 'off');



    % Add a descriptive title
    % titleStr = sprintf('CDF Plot: BW %i MHz, Tx Power 30 dBm, Window %0.2f m', B/1e6, w);
    % title(titleStr, 'Interpreter', 'latex', 'FontSize', 16);

    if length(dim) > 1
        outputFileName = sprintf('%s/FR1_CDF_BW_%iMHz_TxPower30dBm_Window_%0.2fm_3D_fap_threshold_%ddB',output_file_path, B/1e6, w, fap_threshold);
    else
        outputFileName = sprintf('%s/FR1_CDF_BW_%iMHz_TxPower30dBm_Window_%0.2fm_Zaxisfap_threshold_%ddB', output_file_path, B/1e6, w, fap_threshold);
    end
    print(gcf, [outputFileName '.pdf'], '-dpdf', '-r300');  % Save as vector-based PDF
    print(gcf, [outputFileName '.png'], '-dpng', '-r300');  % Save as high-resolution PNG
    savefig(gcf, [outputFileName '.fig']);  % Save as MATLAB figure file
    print(gcf, [outputFileName '.eps'], '-depsc', '-r300');  % Save as EPS
end

%%  FR3
for fap_threshold_idx=1:Nthresh
    fap_threshold = fap_threshold_vec(fap_threshold_idx);

    % Create figure with specific dimensions for publication
    figure('Units', 'centimeters' , 'Position', [1 1 16.8 12], 'PaperPositionMode', 'auto');

    for fileidx=4:6
        % Plot data with specified properties
        legend_str = sprintf("CRLB:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data1(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{1}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx});

        hold on;
        legend_str = sprintf("D-NLS:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data2(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{2}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx}); % Red
        legend_str = sprintf("LLS:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data3(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{5}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx}); % Yellow



        % Customize axes
        xlim([0, edges(end)]);
        ylim([0, 1]);
        grid on;
        xlabel('RMSE (m)', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('CDF', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off');

        ax = gca; % Get the current axes
        ax.XMinorGrid = 'on'; % Enable minor grid for x-axis
        ax.YMinorGrid = 'on'; % Enable minor grid for y-axis

        % Set the spacing for major ticks
        ax.XTick = 0:2:peblimh;
        ax.YTick = 0:0.1:1;

        % Aspect ratio for the plot
        % daspect([5 0.8 1]);
        daspect([10 0.8 1]);

    end
    % Customize the legend for clarity
    legend('Interpreter', 'latex', 'Location', 'northoutside', 'Orientation', 'horizontal','NumColumns', 3,'FontSize', 10,'Box', 'off');

    % Add a descriptive title
    % titleStr = sprintf('CDF Plot: BW %i MHz, Tx Power 30 dBm, Window %0.2f m', B/1e6, w);
    % title(titleStr, 'Interpreter', 'latex', 'FontSize', 16);

    if length(dim) > 1
        outputFileName = sprintf('%s/FR3_CDF_BW_%iMHz_TxPower30dBm_Window_%0.2fm_3D_fap_threshold_%ddB',output_file_path, B/1e6, w, fap_threshold);
    else
        outputFileName = sprintf('%s/FR3_CDF_BW_%iMHz_TxPower30dBm_Window_%0.2fm_Zaxisfap_threshold_%ddB', output_file_path, B/1e6, w, fap_threshold);
    end
    print(gcf, [outputFileName '.pdf'], '-dpdf', '-r300');  % Save as vector-based PDF
    print(gcf, [outputFileName '.png'], '-dpng', '-r300');  % Save as high-resolution PNG
    savefig(gcf, [outputFileName '.fig']);  % Save as MATLAB figure file
    print(gcf, [outputFileName '.eps'], '-depsc', '-r300');  % Save as EPS
end
%% FR2
% Create figure with specific dimensions for publication
for fap_threshold_idx=1:Nthresh
    fap_threshold = fap_threshold_vec(fap_threshold_idx);

    figure('Units', 'centimeters' , 'Position', [1 1 16.8 12], 'PaperPositionMode', 'auto');

    for fileidx=7:9
        % Plot data with specified properties
        legend_str = sprintf("CRLB:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data1(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{1}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx});

        hold on;
        legend_str = sprintf("D-NLS:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data2(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{2}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx}); % Red
        legend_str = sprintf("LLS:%s",frequencies{fileidx});
        plot(edges(1:Nbins), data3(:,fileidx,fap_threshold_idx), 'MarkerSize', 7, 'LineWidth', 2, ...
            'MarkerIndices', markerindices, 'Marker', markervec{5}, ...
            'DisplayName', legend_str, 'Color', colors{fileidx}); % Yellow



        % Customize axes
        xlim([0, edges(end)]);
        ylim([0, 1]);
        grid on;
        xlabel('RMSE (m)', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('CDF', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off');

        ax = gca; % Get the current axes
        ax.XMinorGrid = 'on'; % Enable minor grid for x-axis
        ax.YMinorGrid = 'on'; % Enable minor grid for y-axis

        % Set the spacing for major ticks
        ax.XTick = 0:2:peblimh;
        ax.YTick = 0:0.1:1;

        % Aspect ratio for the plot
        % daspect([5 0.8 1]);
        daspect([10 0.8 1]);

    end
    % Customize the legend for clarity
    legend('Interpreter', 'latex', 'Location', 'northoutside', 'Orientation', 'horizontal','NumColumns', 3,'FontSize', 10,'Box', 'off');

    % Add a descriptive title
    % titleStr = sprintf('CDF Plot: BW %i MHz, Tx Power 30 dBm, Window %0.2f m', B/1e6, w);
    % title(titleStr, 'Interpreter', 'latex', 'FontSize', 16);

    if length(dim) > 1
        outputFileName = sprintf('%s/FR2_CDF_BW_%iMHz_TxPower30dBm_Window_%0.2fm_3D_fap_threshold_%ddB',output_file_path, B/1e6, w, fap_threshold);
    else
        outputFileName = sprintf('%s/FR2_CDF_BW_%iMHz_TxPower30dBm_Window_%0.2fm_Zaxisfap_threshold_%ddB', output_file_path, B/1e6, w, fap_threshold);
    end
    print(gcf, [outputFileName '.pdf'], '-dpdf', '-r300');  % Save as vector-based PDF
    print(gcf, [outputFileName '.png'], '-dpng', '-r300');  % Save as high-resolution PNG
    savefig(gcf, [outputFileName '.fig']);  % Save as MATLAB figure file
    print(gcf, [outputFileName '.eps'], '-depsc', '-r300');  % Save as EPS
end