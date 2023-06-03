function [uHat, stdUHat, lambda2, sigma2Hat, performance] = bayesianSmoothing(sampleData,uTrue,wlen,kstd,scenarioFlag)
%bayesianSmoothing function performs the smoothing of the CGM/raw 
%timeseries, assuming that the measures are corrupted by a non-stationary
%colored noise modelled by an autoregressive model of order 2.
%
% The function handles the presence of missing data (nan) with a virtual
% grid
%
% 
%Input:
%  - sampleData: vector of EGV/raw measurements (raw:pA; egv:mg/dL)
%  - uTrue: vector of ground-truth signal (raw:pA; egv:mg/dL), i.e., IG 
%  (available only in simulated datasets)
%  - wlen: length of the window for the partitioning phase (scalar>20)
%  (samples)
%  - kstd: standard deviation of the gaussian kernel for the reconciliation
%  phase (scalar>0) (samples)
%  - ScenarioFlag: %0=raw data; 1=EGV data
%
%
%Output:
%  - uHat: a vector of the filtered signal (raw: pA; egv:mg/dL)
%  - std: standard deviation of the estimated filtered signal (raw: pA; egv:mg/dL)
%  - lambda2: estimated model error variance (raw: pA^2; egv:mg^2/dL^2)
%  - sigma2Hat: smoothed estimated noise variance (raw: pA^2; egv:mg^2/dL^2)
%  - performance: structure containing the performance metric when u is
%  available
%
% ------------------------------------------------------------------------

% % ------------------------------------------------------------- % %
% % ----------------------- Mirroring --------------------------- % %
% % ------------------------------------------------------------- % %

mirrorHead = sampleData(1:(wlen-1)*2);
mirrorTail = sampleData(end-(wlen-1)*2+1:end);
sampleDataMirrored = [flipud(mirrorHead); sampleData ; flipud(mirrorTail)];


% % ----------------------------------------------------------- % %
% % ---------------------- Initialization --------------------- % %
% % ----------------------------------------------------------- % %

saturationFlag = 1; %1=saturate the estimated lambda2 within a [min,max] range; 0=no saturation
lambda2RangeFlag = 0; %1=population values of [min,max] range; 0 = subject-specific values of [min,max] range;

nSample = length(sampleDataMirrored);
wHalfLen = floor(wlen/2); %Half window length in partitioning phase
kerWeigth = normpdf([-wHalfLen: wHalfLen],0,kstd); %weigths of the gaussian kernel for the reconciliation phase
kerWeigth = kerWeigth/norm(kerWeigth,1); %normalized weights

% Regularization initialization (regularization factor == gamma)
par.gammamin = 1e-10; %Minimum value of regularization factor
par.gammamax = 1e10; %Maximum value of regularization factor
par.max_iter = 1000; %Maximum number of iterations to obtain the regularization parameter
par.tol = 10e-7; %Tolerance



% % ----------------------------------------------------------- % %
% % ---------------- Partitioning & Filtering ----------------- % %
% % ----------------------------------------------------------- % %

[uhat_tmp, cov_u, sigma2, lambda2] = filterPartitionedData(sampleDataMirrored,wlen,par,saturationFlag,lambda2RangeFlag,scenarioFlag);



% % ----------------------------------------------------------- % %
% % -------------- Reconciliation of the signals -------------- % %
% % ----------------------------------------------------------- % %

% Filtered signal reconciliation
j = 1;
for i = wlen : nSample-wlen +1
    uHatSmooth(j) = (diag(fliplr(uhat_tmp(i-wHalfLen:i+wHalfLen,:))))'*kerWeigth';
    stdSmooth(j) = sqrt((diag(fliplr(cov_u(i-wHalfLen:i+wHalfLen,:))))'*kerWeigth');
    j = j+1;
end

% sigma2 reconciliation
j = 1; 
for i = (wHalfLen+1):length(sigma2)-wHalfLen
    sigma2Smooth(j) = ((sigma2(i-wHalfLen : i+wHalfLen))'*kerWeigth');
    j= j+1;
end


% % ------------------------------------------------------------- % %
% % --------------------- De-Mirroring -------------------------  % %
% % ------------------------------------------------------------- % %

uHat = [uHatSmooth(wlen:end-(wlen-1))]';
sigma2Hat = [sigma2Smooth(wlen:end-(wlen-1))];
stdUHat = [stdSmooth(wlen:end-(wlen-1))];


% % ----------------------------------------------------------- % %
% % ----------------- Performance computation ----------------- % %
% % ----------------------------------------------------------- % %
performance = [];
if ~isnan(uTrue)
   % uTrue = uTrue(wlen:end-wlen+1);
    performance.RMSE = rmse(uHat,uTrue); %Root Mean Squared Error
    performance.MARD = mard(uHat,uTrue); %Mean Absolute Relative Difference
end

end