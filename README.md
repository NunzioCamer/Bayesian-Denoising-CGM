# Bayesian-Denoising-CGM
A Bayesian denoising algorithm to deal with colored, non-stationary noise in continuous glucose monitoring timeseries

# Reference
N. Camerlingo, I. Siviero, M. Vettortti, G. Sparacino, S. Del Favero, and A. Facchinetti. "Bayesian denoising algorithm dealing with colored, non-stationary noise in continuous glucose monitoring timeseries", 2023, submitted.

# Description of the algorithm
bayesianSmoothing.m performs smoothing of CGM timeseries, assuming that the measurements are corrupted by non-stationary colored noise, modelled by an autoregressive model of order 2 (Vettoretti et al., Sensors, 2019).

## Input:
- sampleData: CGM measurements (mg/dL)
- uTrue: ground-truth interstitial glucose (mg/dL), available only in simulated scenarios, to compute performance metrics
- wlen: length of the window in the partitioning phase (samples, >20, recommended: 41)
- kstd: standard deviation of the gaussian kernel in the reconciliation phase (samples, >0, recommended: 11)

## Output:
- uHat: denoised glucose signal (mg/dL)
- std: uncertainty around the estimated filtered signal, based on the covariance matrix of the measurment noise (mg/dL)
- lambda2: estimated model error variance (mg^2/dL^2)
- sigma2Hat: estimated noise variance timeseries (mg^2/dL^2)
- performance: structure containing the performance metric when uTrue is available

The algorithm first partitions the CGM timeseries into partially-overlapped windows of a user-defined length (see filterPartitionedData.m).
Then, it performs the Bayesian denoising of each window, i.e., it estimates a regularization parameter which optimally balance the fidelity to data with the roughness of the estimate (De Nicolao et al., Automatica, 1997). This step is achieved using a numerically efficient approach, via singular values decomposition (see filteringAggressivenessEstimation.m). 
Finally, the estimated denoised signals are reconciliated using a smoothing kernel.
The algorithm is robust to missing data (data gaps), handling them using a virtual grid.
