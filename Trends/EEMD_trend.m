% EEMD_trend function copied below
function [TIME, IMFs, R] = EEMD_trend(TIME,TEMP,max_IMF,Trials,noise_std,figon,sift_tol)

%% Remove NaNs

c = isfinite(TEMP);
TEMP = TEMP(c);
TIME = TIME(c);

%% calculate the EEMD

% prepare
I = NaN(numel(TEMP),20,Trials);
R = NaN(numel(TEMP),Trials);

for n = 1:Trials
    % add Gaussian Noise
    noisy_T = (randn(1,numel(TEMP))*noise_std).*nanstd(TEMP) + TEMP;   
    if max_IMF > 0
        [e, R(:,n)] = emd(noisy_T,'MaxNumIMF',max_IMF,'SiftRelativeTolerance',sift_tol);
        I(:,1:size(e,2),n) = e;
    else
        [e, R(:,n)] = emd(noisy_T,'SiftRelativeTolerance',sift_tol);
        I(:,1:size(e,2),n) = e;
    end
end

% Average IMFs and residual
for n = 1:20 
    mn = nanmean(squeeze(I(:,n,:)),2);
    if isfinite(mn) > 0
        IMFs(n,:) = nanmean(squeeze(I(:,n,:)),2);
    end
end
R = nanmean(R,2);


%% plot figure if requested

if figon ==1
    figure('units','normalized','position',[0 0 .5 .9]);
    for n = 1:size(IMFs,1)
        subplot(size(IMFs,1)+1,1,n)
        title(['IMF ', num2str(n)])
        plot(TIME,IMFs(n,:))
        datetick
    end
    subplot(size(IMFs,1)+1,1,n+1)
    title('Trend')
    plot(TIME,R-nanmean(R),'k')
    datetick 
    set(gcf,'Color','W');
    set(gca,'FontSize',10);    
end


end