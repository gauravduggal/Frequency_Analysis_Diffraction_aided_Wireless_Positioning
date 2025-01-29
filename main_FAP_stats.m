clc
clear
close all
root_path = 'Data\';
output_file_path = 'Results\';

%%
if ~exist(output_file_path, 'dir')
    mkdir(output_file_path); % Create the directory
    fprintf('Directory "%s" created.\n', output_file_path);
else
    fprintf('Directory "%s" already exists.\n', output_file_path);
end

%% Loading data


c = physconst('LightSpeed');
base_path_vec = {
    '\Ray Tracing Simulation Code_800MHz_working\0p8GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_2p4GHz_working\2p4GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_5p8GHz_working\5p8GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_2R_2T_2D_APG_2m',
    '\Ray Tracing Simulation Code_8GHz_working\8GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_10GHz_working\10GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_15GHz_working\15GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m',
    '\Ray Tracing Simulation Code_28GHz_working\28GHz_4TX_30dBm_Combined_Data_4Anchors_diffused_mpc_2R_2T_2D',
    '\Ray Tracing Simulation Code_37GHz_working\37GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m'
    '\Ray Tracing Simulation Code_48GHz_working\48GHz_30dBm_Combined_Data_4Anchors_diffused_mpc_6R_6T_1D_APG_5m'
    };

Nf = length(base_path_vec);
%number of anchors
Na = 4;
base_path = base_path_vec{1};
load(sprintf('%s%s%s',root_path,base_path,'\Transmitter1_Combined.mat'),'data');
tx_data(:,1) = data;

%System parameters
Nr = size(tx_data,1);
fap_threshold_vec = [0,10,20];
B = 400e6;
k = physconst('Boltzmann');
T = 300;
%Noise Floor
N0 = k*T*B;

%FAP statistics
FAP_SNR = zeros(Nr,Na,Nf);
rx_processing_gain_vec=[-10,-10,-10,0,0,0,20,20,20];
fap_class = zeros(Nf,4);
for fap_threshold_idx = 1:length(fap_threshold_vec)

    fap_threshold = fap_threshold_vec(fap_threshold_idx);

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

        clear data





        edges = 80:1:110;
        clear data
        Nr = length(tx_data(:,1));
        rx_loc = get_receiver_locations(tx_data(:,1));

        % Anchor positions
        tx_loc(:,1) = get_transmitter_location(tx_data(:,1));
        tx_loc(:,2) = get_transmitter_location(tx_data(:,2));
        tx_loc(:,3) = get_transmitter_location(tx_data(:,3));
        if Na==4
            tx_loc(:,4) = get_transmitter_location(tx_data(:,4));
        end

        %anchor positions
        a = zeros(size(tx_loc));
        %anchor coordinates transformed from WI to TWC paper model
        %y coordinate (LHS)
        a(2,:) = tx_loc(1,:);
        %x coordinate (LHS)
        a(1,:) = tx_loc(2,:);
        %z coordinate (LHS)
        a(3,:) = tx_loc(3,:);


        % Pe for Tx-D-Rx, Tx-D-T-Rx, Tx-D-T-T-Rx
        Pe_mpc3 = 0;
        % Pe for Tx-Rx,
        Pe_mpc1 = 0;
        % Pe for Tx-R-Rx,
        Pe_mpc2 = 0;
        mpc3 = {'Tx-D-Rx','Tx-D-T-Rx','Tx-D-T-T-Rx','Tx-D-T-T-T-Rx','Tx-D-T-T-T-T-Rx'};
        mpc1 = {'Tx-Rx','Tx-T-Rx','Tx-T-T-Rx','Tx-T-T-T-Rx','Tx-T-T-T-T-Rx','Tx-T-T-T-T-T-Rx',};
        mpc2 = {'Tx-R-Rx','Tx-R-T-Rx','Tx-R-T-T-Rx','Tx-R-T-T-T-Rx','Tx-R-T-T-T-T-Rx','Tx-R-T-T-T-T-T-Rx'};
        temp_fap_class = [0,0,0,0];
        for ridx = 1:Nr
            %node position
            xn = rx_loc(2,ridx);
            yn = rx_loc(1,ridx);
            zn = rx_loc(3,ridx);
            np = [xn;yn;zn];
            for aidx=1:Na
                if tx_data(ridx,aidx).NumberOfPaths==0
                    continue
                end
                %FAP statistics
                %pick first arriving path for rx location ridx
                FAP = get_first_receiver_path(tx_data(ridx,aidx),fap_threshold);
                %calculated signal power in FAP
                FAP_rxPowerdBm = FAP.ReceivedPower_dBm_ + rx_processing_gain;
                FAP_rxPower = 10^((FAP_rxPowerdBm)/10)*1e-3;
                FAP_snr_linear = FAP_rxPower/N0;
                FAP_snr_db = 10*log10(FAP_snr_linear);
                FAP_SNR(ridx,aidx, base_path_idx) =  FAP_snr_db;
             
                temp_fap_class = classify_FAP_mpc_group(FAP, mpc1,mpc2,mpc3,temp_fap_class);
            end
        end
       


        fap_class(base_path_idx,:) = 100*temp_fap_class/(Na*Nr);

    end
   
    %% P_FAP Statistics
    % Set up the figure with dimensions optimized for publication
    frequencies = {'0.8', '2.4','5.8','8','10','15', '28','37','48'};
    figure('Units', 'inches', 'Position', [1 1 6 3.75], 'PaperPositionMode', 'auto');



    % Create the stacked bar plot
    bar(fap_class, 'stacked','BarWidth', 0.5);

    xline(3.5, '--k', 'LineWidth', 2.5);
    xline(6.5, '--k', 'LineWidth', 2.5);
    text(1.75, 105, '$FR1$', 'HorizontalAlignment', 'center', FontSize=14, Interpreter='latex');
    text(4.75, 105, '$FR3$', 'HorizontalAlignment', 'center', FontSize=14,Interpreter='latex');
    text(7.75, 105, '$FR2$', 'HorizontalAlignment', 'center', FontSize=14,Interpreter='latex');
    set(gca, 'XTickLabel', frequencies, 'FontSize', 12, 'FontWeight', 'bold'); % Add frequency labels

    % Customize the axes for readability and research standards
    ylabel('$P_{FAP} (\%)$', 'FontSize', 14, 'Interpreter', 'latex');

    xlabel('Frequency (GHz)', 'FontSize', 14, 'Interpreter', 'latex');

    ylim([0, 110]); % Set the y-axis range to 0-100%
    grid on; % Enable grid lines
    yticks(0:10:100); % Set y-ticks at intervals of 10
    % yticklabels(arrayfun(@num2str, 0:10:100, 'UniformOutput', false));

    % Use LaTeX interpreter for the y-axis label

    % Add a legend for clarity, ensuring font size and position are suitable
    legend({'MPC-1', 'MPC-2', 'MPC-3', 'MPC-4'}, ...
        'Location', 'northoutside', 'Orientation', 'horizontal', ...
        'FontSize', 12, 'Box', 'off', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');

    % Add minor grid lines for improved readability (optional)
    set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'on');

    %Set FR1, FR2, FR3 regions


    % Adjust the figure background and other aesthetics
    set(gcf, 'Color', 'w'); % Set figure background to white for publication

    % Generate a dynamic file name for saving
    outputFileName = sprintf("%s/FAP_stat_BW_%d_MHz_fap_threshold_%ddB", ...
        output_file_path, B / 1e6, fap_threshold_vec(fap_threshold_idx));

    % Save the figure in high-quality formats for publications
    print(gcf, sprintf("%s.pdf", outputFileName), '-dpdf', '-r300');  % Save as vector-based PDF
    print(gcf, sprintf("%s.png", outputFileName), '-dpng', '-r300');  % Save as high-resolution PNG
    print(gcf, sprintf("%s.eps", outputFileName), '-depsc', '-r300');  % Save as EPS

    %% FAP SNR boxplot
    Nf = size(FAP_SNR,3);
    bands = {'FR1','FR3','FR2'};

    figure('Units', 'inches', 'Position', [1 1 6 3.75], 'PaperPositionMode', 'auto');
    % Define frequency labels for the x-axis

    % Create a high-resolution figure for publication

    %FAP SNR boxplot
    for idx = 1:Nf

        data = FAP_SNR(:,idx);

        % Calculate specific percentiles
        p5 = prctile(data, 4);   % 5th percentile
        p25 = prctile(data, 25); % 25th percentile (Q1)
        p50 = prctile(data, 50); % 50th percentile (median)
        p75 = prctile(data, 75); % 75th percentile (Q3)
        p90 = prctile(data, 90); % 90th percentile


        % Plot the interquartile range (box)
        fill([idx-0.1 idx+0.1 idx+0.1 idx-0.1], [p25 p25 p75 p75], [0 0.4470 0.7410], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);
        hold on;

        % Plot the median line
        plot([idx-0.1 idx+0.1], [p50 p50], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5);

        % Plot whiskers
        plot([idx idx], [p5 p25], 'k--', 'LineWidth', 1.2); % Lower whisker
        plot([idx idx], [p75 p90], 'k--', 'LineWidth', 1.2); % Upper whisker
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gca, 'XTickLabel', frequencies, 'FontSize', 12, 'FontWeight', 'bold');
        % Customize labels and title
        ylabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', 14);
        xlabel('Frequency (GHz)','Interpreter', 'latex', 'FontSize', 14);
        % title('SNR Distribution Analysis', 'Interpreter', 'latex', 'FontSize', 16);
        % set(gca, 'FontSize', 12, 'Box', 'off', 'LineWidth', 1.2,'Interpreter', 'latex');
        ylim([-40,60])

        % Add minor grid lines for improved readability (optional)
        set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'on');
        % Set axis limits and adjust aspect ratio
        % xlim([0.5 1.5]);
        % ylim([p5, p90]);
        % aspect = 5;
        % tempy = p90 - p5;
        % tempx = 1;
        % daspect([tempx tempy/aspect 1]);
        yticks(-40:10:60);

        % Add grid for better readability
        grid on;
    end
    xline(3.5, '--k', 'LineWidth', 2.5);
    xline(6.5, '--k', 'LineWidth', 2.5);
    text(1.75, 55, '$FR1$', 'HorizontalAlignment', 'center', FontSize=14, Interpreter='latex');
    text(4.75, 55, '$FR3$', 'HorizontalAlignment', 'center', FontSize=14,Interpreter='latex');
    text(7.75, 55, '$FR2$', 'HorizontalAlignment', 'center', FontSize=14,Interpreter='latex');
    set(gca, 'XTickLabel', frequencies, 'FontSize', 12, 'FontWeight', 'bold'); % Add frequency labels

    xticks(1:(Nf+1))
    xlim([0,(Nf+1)])
    xticklabels(frequencies)
    outputFileName = sprintf("%s\\FAP_SNR_BW_%d_MHz_fap_threshold_%ddB",output_file_path,B/1e6,fap_threshold_vec(fap_threshold_idx));
    print(gcf, sprintf("%s.pdf",outputFileName), '-dpdf', '-r300');  % Save as vector-based PDF
    print(gcf, sprintf("%s.png",outputFileName), '-dpng', '-r300');  % Save as high-resolution PNG
    print(gcf, sprintf("%s.eps",outputFileName), '-depsc', '-r300');  % Save as EPS


end