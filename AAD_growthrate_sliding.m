function dlogOD = AAD_growthrate_sliding(time,logOD,step)

dimtime = find(size(logOD)==length(time));
if isempty(dimtime)
    error('dimension of time vector does not match rows or columns in logOD')
elseif length(dimtime)==2
    warning('relabs is square, assuming that rows are time columns are samples')
    dimtime = 1;
end

if dimtime==2
    logOD = logOD';
end

dlogOD = nan(size(logOD));
for i = 1:size(logOD,2)
    for k = 1:length(time)
        x = [ones(length(max(k-floor(step/2),1): min(k+floor(step/2),length(time))),1),time(max(k-floor(step/2),1):min(k+floor(step/2),length(time)))];
        y = logOD(max(k-floor(step/2),1):min(k+floor(step/2),length(time)),i);
        temp = x\y;
        dlogOD(k,i) = temp(2);
    end
end