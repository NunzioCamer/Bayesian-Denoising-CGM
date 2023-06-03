function [uHat,covUHat,sigma2,lambda2] = filteringAggressivenessEstimation(sampleData, noiseCorrPar, par, regularizationCriterion, fixedLambda2)
% filteringAggressivenessEstimation implements the regularization criteria 
% using singular values decomposition (SVD) assuming that the measures are 
% corrupted by a non-stationary colored noise modelled by an autoregressive 
% model of order 2. The regularization criteria are the so-called 
% "Consistency criterion #3" if sigma2 and lambda2 are unknown, or the 
% "Consistency criterion #2" if sigma2 is unknown but lambda2 is fixed.
%
% The function handles the presence of missing data (nan) using a virtual
% grid
% 
% Input:
%  - sampleData: vector of EGV/raw measurements (raw:pA; egv:mg/dL)
%  - noiseCorrPar: vector of parameters of the autocorrelation model
%  describing the measurement noise (Possible value: [-1.3, 0.42] for AR(2))
%  - par: hyperparameters for the regularization pipeline to find the
%  optimal filtering aggressiveness
%  - regularizationCriterion: criterion for the estimation of the filtering
%  aggressiveness: 2=consistency criterion #2 (i.e., sigma2 to be
%  estimated, lambda2 fixed); 3=consistency criterion #3 (i.e., sigma2 to
%  be estimated, lambda2 to be estimated) 
%  - lambda2RangeFlag: %1=population values of [min,max] range;
%  0 = subject-specific values of [min,max] range;
%  - ScenarioFlag: %0=raw data; 1=EGV data
%
% Output:
%  - uHat: a vector of the filtered signal (raw:pA; egv:mg/dL)
%  - covUHat: covariance matrix of the estimated filtered signal (raw:pA^2;
%  egv:mg^2/dL^2)
%  - sigma2: estimated value of variance of the measurement noise (raw:pA^2; 
%  egv:mg^2/dL^2)
%  - lambda2: estimated value of variance of the model error (raw:pA^2;
%  egv:mg^2/dL^2)
%
% 
% References:
% De Nicolao et al., "Non parametric input estimation in a physiological 
%           system: problems, methods and case studies", Automatica, 1997
% Sparacino et al., "A Bayesian approach to estimate evoked potentials",
%           Comput. Methods Programs. Biomed., 2002
% 
% ------------------------------------------------------------------------

% % ----------------------------------------------------------- % %
% % ---------------- Virtual grid definition ------------------ % %
% % ----------------------------------------------------------- % %

[idxNan,~] = find(isnan(sampleData));     %find missing samples
N = length(sampleData);
sampleData(idxNan) = [];                  %timeseries without nan 
n = length(sampleData);

vt = [1:1:N];
vt(idxNan)=[];

Gv = eye(N);
G = Gv(vt',:);


% % ----------------------------------------------------------- % %
% % ---------------- Parameters inizialization ---------------- % %
% % ----------------------------------------------------------- % %

%SMOOTHING
a1 = noiseCorrPar(1); %noise correlation parameters
a2 = noiseCorrPar(2);

% Initialization of A (describes the correlation of measurement error)
A_col = zeros(n,1);
A_col(1) = 1;
A_col(2) = a1;
A_col(3) = a2;

A_row(1) = 1;
A_row(1,n) = 0;
A = toeplitz(A_col,A_row);

% matrix F
m = 2; %2 for double-integrated random walk, to model IG
F_col = zeros(N,1);
F_col(1) = 1;
invF_col = F_col;
for k=1:m
    F_col = filter([1 -1],1,F_col);
    invF_col = filter(ones(N,1),1,invF_col);
end
F_row(1) = 1;
F_row(1,N) = 0;
F = toeplitz(F_col,F_row); %Measures the roughness of the estimate



% % ----------------------------------------------------------- % %
% % -------------- Singular Values Decomposition -------------- % %
% % ----------------------------------------------------------- % %
B = inv(A'*A);
H = B^(-1/2)*G*inv(F);
[U,D,V] = svd(H);
diagD = diag(D);
DtraspD = diagD.^2;

xi = U'*B^(-1/2)*(sampleData-sampleData(1));


% % ----------------------------------------------------------- % %
% % ----------------- Regularization Criteria ----------------- % %
% % ----------------------------------------------------------- % %
gammamin = par.gammamin;
gammamax = par.gammamax;
max_iter = par.max_iter;
toll = par.tol;

itergamma= 1; %number of iterations
reachedConvergence=0;


while reachedConvergence==0
   % update gamma
    gamma = 10^((log10(gammamax)+log10(gammamin))/2);
    % estimation in the new coordinates
    eta = diagD./(DtraspD+gamma).*xi;
    eta(n+1:N)=0;
    % residuals in the new coordinates
    ro = gamma*xi./(DtraspD+gamma);
    % degrees of freedom in the new coordinates
    qgamma = sum(diagD.^2./(DtraspD+gamma));
    wrss = sum(ro.^2);
    wess = sum((diagD.*xi./(DtraspD+gamma)).^2);

  switch regularizationCriterion 
      case 3
        wrss_term = wrss/(N-qgamma);
        wess_term = gamma*wess/qgamma;
        difference = (wrss_term-wess_term); %Consistency criterion #3
        diff = abs(difference);
      case 2
        sigma2 = fixedLambda2*gamma;
        difference = wrss-sigma2*(N-qgamma); %Consistency criterion #2
        diff = abs(difference/wrss); 
  end

    if difference<0
        gammamax=gamma; % reduce regularity
    else
        gammamin=gamma; % increase regularity
    end %if

    if diff<toll || itergamma==max_iter
        reachedConvergence=1; 
    end %if
    itergamma=itergamma+1;
  
end

eta = (diagD.*xi)./(DtraspD+gamma);
eta(n+1:N) = 0;
qgamma = sum(DtraspD./(DtraspD+gamma));

C = filter(invF_col,1,V);

uHat = C*eta;
uHat = uHat+sampleData(1);  %Estimated signal 

switch regularizationCriterion 
      case 3 %consistency criterion #3
        sigma2=wrss/(N-qgamma); %estimated noise variance
        lambda2=sigma2/gamma;   %estimated model error variance
      case 2 %consistency criterion #2
        sigma2=fixedLambda2*gamma; %estimated noise variance
        lambda2 = fixedLambda2;    %fixed model error variance
end

covEta = sigma2*ones(n,1)./(DtraspD+gamma);
covEta(n+1:N,1)=lambda2;

c2 = C.^2;
covUHat = c2*covEta;  %Estimated covariance

end