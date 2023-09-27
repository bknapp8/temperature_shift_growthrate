%% Test script for extracting response time using a 'within error' constraint
close all;
figure; 
delta_gr = gr - gr_final;
plot(time_shifted, delta_gr, 'k')
xlabel('Time (min)')
ylabel('G - G_{ss}')
mintime = 15; % Minimum time to consider for response time

% Plot error lines around zero point
hold on;
plot(time_shifted,repmat(-gr_final_sem, length(time_shifted),1), 'r')
plot(time_shifted,repmat(gr_final_sem, length(time_shifted),1), 'r')

% Extract response time
ResponseTime = min(time_shifted(abs(delta_gr)<0.005)); % Time at which growth rate is within 0.5% of steady state
ResponseTime_err = abs(ResponseTime - min(time_shifted(abs(delta_gr)<gr_final_sem)));