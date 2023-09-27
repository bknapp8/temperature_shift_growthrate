figure;
plot(time_norm,gr_norm);
tcut_min = 0.25;
tcut_max = min(time_norm(gr_norm>0.95));
tfit = time_norm(time_norm<tcut_max & time_norm>tcut_min);
gfit = gr_norm(time_norm<tcut_max & time_norm>tcut_min);
efit = err_norm(time_norm<tcut_max & time_norm>tcut_min);
hold on;
plot(tfit,gfit, 'k','LineWidth',1);
box off;

%% Perform weighted fitting
f = fitlm(tfit,gfit,'Weights',efit);
xeval = linspace(0,3,100);
hold on;
plot(xeval,f.feval(xeval));
xlim([-1 3]);
ylim([0 1.5])

%% Now evaluate response time, with error
yint = f.Coefficients.Estimate(1);
yint_err = f.Coefficients.SE(1);
slope = f.Coefficients.Estimate(2);
slope_err = f.Coefficients.SE(2);

RT = (1 - yint)/slope;
RT_err = sqrt((yint_err/slope)^2 + RT^2*(slope_err/slope)^2);
RT_abs = RT*log(2)/gr_final*60;
RT_abs_err = RT_err*log(2)/gr_final*60;
