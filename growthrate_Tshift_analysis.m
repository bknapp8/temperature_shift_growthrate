%% Specify T-shift metadata and options
close all;
time_shift = 13.5; % Time at which T-shift occurs (unit = minutes)
normalizing = true; % Should normalization analysis take
saving = true;
prefix = '20191208_BK001_LB_Tshift_33to37'; % Prefix for saving figures
filename = '20220323_BK001_LB_Tshift_27to37_feature_0_24-Mar-2022_CONTOURS_TOTALS.mat';
load(filename);

%% Define parameters for filtering
traj_total = TOTALS; % Structure containing all trajectories
t_int = 0.5; % Frame interval (min.)
res = 0.11; % Pixel resolution (microns)
min_traj_len = 10; % Minimum trajectory length
begcut = 4; % Initial frames of trajectory to remove
endcut = 0; % Final frames of trajectory to remove
smf = 4;
Rmin = -.1;
Rmax = 3;
plot_toggle = 1;

%% Perform filtering
[traj_filtered, T_filtered, R_filtered, Ncell] = filter_trajectories(traj_total, ...
    t_int, res, min_traj_len, begcut, endcut, smf, Rmin, Rmax, plot_toggle);

%% Plot binned rates using Standard Error of the Mean

% Bin data
gx = (min(T_filtered):t_int*1:max(T_filtered))';
[b,n,s] = bindata(T_filtered,R_filtered, gx);

% Plot shaded error bar
f1 = figure;
shadedErrorBar(gx,b,s./sqrt(n),'lineProps',{'-k', 'linewidth', 2},'transparent',1);
set(gca, 'fontsize', 20);
xlabel('Time (min)');
ylabel('Growth rate (1/h)');
set(gcf, 'position', [0 0 400 300])
box off;
if saving
    savefig(f1, strcat(prefix,'_growthrate_sm',num2str(smf),'_N', num2str(Ncell), '_SEM'));
end

%% Plot Shaded Error (SEM) with T-shift correction (time zero = T-shift)

% Plot shaded error bar
f2 = figure;
shadedErrorBar(gx-time_shift,b,s./sqrt(n),'lineProps',{'-k', 'linewidth', 2},'transparent',1);
set(gca, 'fontsize', 20);
xlabel('Time after shift (min)');
ylabel('Growth rate (1/h)');
set(gcf, 'position', [0 0 400 300])
box off;
if saving
    savefig(f2, strcat(prefix,'_growthrate_sm',num2str(smf),'_N', num2str(Ncell), '_SEM_shift'));
end


%% Plot normalized data (if normalizing) with SEM error bars
gx_0 = gx-time_shift;
frames_at_steadystate = 30; % Number of frames at steady state (end of trajectory)

time_shifted = gx - time_shift;
gr = b;
gr_final = mean(b(end-frames_at_steadystate:end)); % NOTE: This sets the final, steady-state value for the growth rate
gr_final_std = sum(s(end-frames_at_steadystate:end).^2)/frames_at_steadystate;
gr_final_sem = gr_final_std./sqrt(frames_at_steadystate);
tau_double = log(2)/gr_final*60;
gr_initial = b(gx_0==0);
gr_norm = (b-gr_initial)/(gr_final-gr_initial);
err_norm = s/(gr_final - gr_initial);
time_norm = gx_0/tau_double;

% Plot normalized data with max growth rate 
f3 = figure;
% Label maximum growth rate in graph
time_test = linspace(-1,4,100);
gr_test = repmat(1, 100,1);

shadedErrorBar(time_norm,gr_norm,err_norm./sqrt(n),'lineProps',{'-k', 'linewidth', 2},'transparent',1);
hold on;
plot(time_test, gr_test, 'k');
xlim([-1 3]);
set(gcf, 'Position', [0 0 400 300]);
set(gca, 'fontsize', 20);
box off;
xlabel('Normalized time (doubling at T final)');
ylabel('Normalized growth rate (rate at T final)');
if saving
    savefig(f3, strcat(prefix,'_growthrate_sm',num2str(smf),'_N', num2str(Ncell), '_SEM_normalized0to1'));
end

