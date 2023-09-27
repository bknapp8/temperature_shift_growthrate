function [traj_filtered, T_filtered, R_filtered, Ncell] = filter_trajectories(traj_total, t_int, res, min_traj_len, begcut, endcut, smf, Rmin, Rmax, plot_toggle);
%% Input parameters
% 1. traj_total = output from cell_tracking_bacteria analysis (usually called
% TOTALS)
% 2. t_int = time interval between frames (unit = minutes)
% 3. res = image resolution (microns/pixel)
% 4. min_traj_len = minimum trajectory length (recommended = 10)
% 5. begcut = cut-off for removing initial frames (recommended = 4)
% 6. endcut = cut-off for removing final frames of trajectory (recommended
% = 0)
% 7. smf = smoothing factor for filtering trajectories (recommended = 4)
% 8. Rmin = minimum growth rate allowed in filtered trajectory (recommended = -0.1)
% 9. Rmax = maximum growth rate allowed in filtered trajectory (this may
% require some checking, but recommended is 1.5*average
% 10. plot_toggle = binary for choosing whether to plot data (0 or 1)

%% Output parameters
% 1. traj_filtered = structure containing filtered trajectories
% 2. T_filtered = total vector of filtered time points (unit = minutes)
% 3. R_filtered = total vector of filtered growth rates (unit = 1/hour)

tree = traj_total;
TOTALS_filtered = [];

T_filtered = [];
R_filtered = [];
figure;

ind_keep = [];
Ncell = 0;
for k=1:length(tree)
    clearvars L T;
    if length(tree(k).traj)>min_traj_len
        for j=1:length(tree(k).traj)
            L(j) = tree(k).traj(j).length;
            T(j) = tree(k).traj(j).frame_num;
        end

        % Smooth length data and obtain growth rate
        Lsm = smooth(smooth(L(begcut:end-endcut),smf),smf)*res;
        R = AAD_growthrate_sliding(T(begcut:end-endcut)'*t_int, log(Lsm), smf)*60; 
        
        % Perform quality cuts on data
        Rmin = -.1;
        Rmax = 3;
        Rdiffmax = 1;
        if max(R)<Rmax & min(R)>Rmin & sum(abs(diff(R))>Rdiffmax)<1
            ind_keep = [ind_keep; k];
            Ncell = Ncell+1;

            T_filtered = [T_filtered; (T(begcut:end-endcut))'*t_int];
            R_filtered = [R_filtered; R];
            
            % Plot if chosen
            if plot_toggle
                subplot(1,2,1);
                plot(T(begcut:end-endcut)*t_int, log(Lsm), 'r');
                if false
                    text(T(end-endcut)*t_int, log(Lsm(end)), num2str(k));
                end
                set(gca, 'fontsize', 20);
                xlabel('Time (min)');
                ylabel('Log(length)');
                title('Filtered cell lengths')
                
                hold on;
                subplot(1,2,2);   
    
                plot(T(begcut:end-endcut)*t_int, R,'r');
                hold on;
                
                set(gca, 'fontsize', 20);
                xlabel('Time (min)');
                ylabel('Growth rate (1/h)');
                title('Filtered growth rates')
                set(gcf, "Position", [0 0 800 300]);
            end
            
        end
    end
end
disp('Number of trajectories analyzed:');
disp(Ncell);

traj_filtered = traj_total(ind_keep);
end
