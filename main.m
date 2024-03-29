%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------%
%            Task 4.1 - Noise characterization and reduction             %
%------------------------------------------------------------------------%
%                    Provided by UNIPD, Sep. 2021                        %
%------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear; 
close all; 
clc

% % ----------------------------------------------------------- % %
% % ----------------------- Load data ------------------------- % %
% % ----------------------------------------------------------- % %

addpath(genpath('src'))
CGMdata = load ('Camerlingo_simulated_2023.mat'); 
CGMdata = CGMdata.simulatedData;
subj = 1;  % selected trace
id_subj = strrep(string(CGMdata(subj).UniqueID),'_',''); % ID subject

firstNotNan = find(~isnan(CGMdata(subj).EGV.glucose),1,'first'); %first sample after calibration
lastNotNan = find(~isnan(CGMdata(subj).EGV.glucose),1,'last');

sampleData = CGMdata(subj).EGV.glucose(firstNotNan:lastNotNan);

time = datetime(datestr(CGMdata(subj).EGV.TimeStamp(firstNotNan:lastNotNan)));  %TimeStamp
uTrue = nan; %ground-truth signal


% % ----------------------------------------------------------- % %
% % ----------------- Hyperparameters setting  ---------------- % %
% % ----------------------------------------------------------- % %

wlen = 41; %length of CGM window for the partitioning phase (samples)
kernel_std = 10; %standard deviation of the gaussian kernel in the 
                %reconciliation phase (samples)


% % ------------------------------------------------------------ % %
% % -------------------- Bayesian smoothing -------------------- % %
% % ------------------------------------------------------------ % %
[uHat, stdUHat, lambda2, sigma2] =  bayesianSmoothing(sampleData,uTrue,wlen,kernel_std,1);

% % ------------------------------------------------------------- % %
% % ----------------------- Plot results  ----------------------- % %
% % ------------------------------------------------------------- % %

figure
ax(1) = subplot(3,1,1:2);
plot(time,sampleData,'k.-','linewidth',1); hold on;
plot(time,uHat,'r','LineWidth',1.5)

fill([time; flipud(time)],[(uHat)+stdUHat' ; flipud((uHat)-stdUHat')],'y','FaceAlpha',0.3);
legend('Sample data','Filtered data','CI (+/-SD)')
axis tight;

ax(2) = subplot(3,1,3);
plot(time, sigma2','r','LineWidth',2)
title('Estimated noise variance')
xlabel('[MM/DD/YY]')

linkaxes(ax,'x')

subplot(3,1,1:2)
ylabel('[mg/dL]')
title(['Subject ', char(id_subj),': EGV data'])
set(gca,'Fontsize',12);

subplot(3,1,3)
ylabel('[mg^2/dL^2]')
set(gca,'Fontsize',12);
