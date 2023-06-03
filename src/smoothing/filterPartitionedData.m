function [uHat, covUHat, sigma2, lambda2] = filterPartitionedData(sampleData,wlen,par,saturationFlag,lambda2RangeFlag,scenarioFlag)
% filterPartitionedData function first partions the signal into windows 
% (fixed user-defined length equal to "wlen") and then filters the data
% (with/without lambda saturation)
%
% The function handles the presence of missing data (nan) using a virtual
% grid
% 
% Input:
%  - sampleData: vector of EGV/raw measurements (raw:pA; egv:mg/dL)
%  - wlen: length of the window for the partitioning phase (scalar>20)
%  (samples)
%  - par: hyperparameters for the regularization pipeline to find the
%  optimal filtering aggressiveness
%  - saturationFlag: %1=saturate the estimated lambda2 within a [min,max] 
%  range; %0=no saturation
%  - lambda2RangeFlag: %1=population values of [min,max] range;
%  0 = subject-specific values of [min,max] range;
%  - ScenarioFlag: %0=raw data; 1=EGV data
%
%
% Output:
%  - uHat: a matrix containing the estimated signals for each window 
%  [n_sample x wlen]
%  - covUHat: covariance matrix of the estimated filtered signal (raw:pA^2;
%  egv:mg^2/dL^2) [n_sample x wlen].
%  - sigma2: vector of the estimated values of variance of the measurement
%    error (raw:pA^2; egv:mg^2/dL^2) [(n_sample-2*(wind_len/2)) x 1]
%  - lambda2: vector of the estimated values of variance of the model 
%  error (raw:pA^2; egv:mg^2/dL^2) [(n_sample-2*(wlen/2)) x 1]
% 
% ------------------------------------------------------------------------


% % ----------------------------------------------------------- % %
% % ---------------------- Initialization --------------------- % %
% % ----------------------------------------------------------- % %

noiseCorrPar = [-1.3, 0.42]; %parameters of the AR(2) describing the measurement noise [Vettoretti et al., Sensors, 2019] 
nSample = length(sampleData);
wHalfLen = floor(wlen/2); %Half window length in partitioning phase
sigma2 = [];
lambda2 = [];
uHat = nan(nSample,wlen);
covUHat = nan(nSample,wlen);


% % ----------------------------------------------------------- % %
% % -------------- Filtering partitioned windows -------------- % %
% % ----------------------------------------------------------- % %

for idxWindStart = (wHalfLen+1) : nSample-wHalfLen 

    yWind = sampleData(idxWindStart-wHalfLen:idxWindStart+wHalfLen); %data partition
    n_gap = length(find(isnan(yWind)));
    
    if n_gap >= wlen-1
        uHat(idxWindStart,:) = nan;
        sigma2Hat = nan;
        lambda2Hat = nan;
    
    elseif n_gap > ceil(0.1*wlen) 
        switch scenarioFlag
            case 1 %EGV data
                lambda2Mean= 1.3047;
                
            case 0 %raw data
                lambda2Mean= 1433.6;
        end
        [uHat(idxWindStart,:),covUHat(idxWindStart,:),sigma2Hat,lambda2Hat] = filteringAggressivenessEstimation(yWind,noiseCorrPar,par,2,lambda2Mean);  %consistency 2
    else
        
        [uHat(idxWindStart,:),covUHat(idxWindStart,:),sigma2Hat,lambda2Hat] = filteringAggressivenessEstimation(yWind,noiseCorrPar,par,3);  %consistency 3
    end
       
    sigma2 = [sigma2; sigma2Hat];
    lambda2 = [lambda2; lambda2Hat];

end


% % ----------------------------------------------------------- % %
% % -------------------- Lambda2 saturation ------------------- % %
% % ----------------------------------------------------------- % %

switch saturationFlag
    case 0 %without saturation of estimated lambda2
        % DO ANYTHING

    case 1 %saturates the values of lambda2 in a [min,max] range
        
        switch lambda2RangeFlag
            case 0 %lambda2 [min,max] range selected as subject-specific values
                
                lambda2Min = prctile(lambda2,10);   %individual values 
                lambda2Max = prctile(lambda2,90);   %individual values
                
            case 1 %lambda2 [min,max] range selected as population values
                switch scenarioFlag
                    case 1 %EGV data
                        lambda2Min = 0.004;
                        lambda2Max = 14.9268;
                    case 0 %raw data 
                        lambda2Min =  17.2586;
                        lambda2Max = 9303.2;
                end
        end
        sigma2 = [];
        lambda2 = [];

        for idxWindStart = (wHalfLen+1 ) : nSample-wHalfLen 
            yWind = sampleData(idxWindStart-wHalfLen:idxWindStart+wHalfLen); %data partition

            % Filtering aggressiveness estimation
            [uHat(idxWindStart,:),covUHat(idxWindStart,:),sigma2Hat,lambda2Hat] = filteringAggressivenessEstimation(yWind,noiseCorrPar,par,3);  %consistency 3

            if lambda2Hat > lambda2Max
                lambda2Sat = lambda2Max;
                [uHat(idxWindStart,:),covUHat(idxWindStart,:),sigma2Hat,~] = filteringAggressivenessEstimation(yWind,noiseCorrPar,par,2,lambda2Sat); %consistency 2
                lambda2Hat = lambda2Sat;
            elseif lambda2Hat < lambda2Min
                lambda2Sat = lambda2Min;
                [uHat(idxWindStart,:),covUHat(idxWindStart,:),sigma2Hat,~] = filteringAggressivenessEstimation(yWind,noiseCorrPar,par,2,lambda2Sat); %consistency 2
                lambda2Hat = lambda2Sat;
            end

            sigma2 = [sigma2; sigma2Hat];
            lambda2 = [lambda2; lambda2Hat];
        end


end