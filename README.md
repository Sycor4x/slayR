# R-convenience-functions
Several functions that are useful for statistical analysis, machine learning and optimization.

**slayR: Bayesian Optimization**

At present, this repo only implements a global minimization algorithm via the function slayR. It implements the Bayesian Optimization algorithm with full hyperparameter marginalization via Hamiltonian Monte Carlo. The package rstan is used for that. slayR implements two acquisition functions, expected improvement (EI) and lower confidence bound (LCB). Optimization of the acquisiton function is done at random via latin hypercube sampling from the posterior. This might seem ad hoc but experiments provide evidence that this approach is only slightly less effective than more expensive methods.

**ML integration**

Eventually, the goal is to write a series of functions which will make the nested cross-validation of large machine learning models easier and quicker. This will be accomplished by writing wrapper methods integrating ML models -- SVMs, neural networks, etc -- into the slayR method. At present, only the SVM component of this is fleshed out to any degree.

Goals for future releases of slayR: 
  - additional acquisition functions (specificaly probability of improvement and at least one entropy-based methods) as well as portfolio acquisition and predictive entropy search
  - more sophisticated methods for maximizing the acquisition function
  - automated input/output warping (presently this is entirely left to the user)
  
Goals for future releases of ML functions integrated into slayR:
  - complete SVM wrapper functions
  - implement wrapper functions and appropriate GP kernel functions for ANNs and regularized regresison models
  
Goals for regularized regression models:
  - implement a regularization model that uses hyperparameter marginalization & compare results to GLMNET

**Spatial Prediction Tasks**

Some of the tasks we work on require spatial interpolation -- for example, some sensors are fixed at specific locations, and we'd like to have "best guesses" of the measurements of the sensor readings at points between those locations. We can use very similar technology to the global optimization setting to create high-quality predictions of the sensor readings at those locations. The methods presently implemented are premised on solely having location data and sensor readings in hand, but simulation experiments demonstrate that this is sufficient to reconstruct sensor readings in a variety of conditions.

As implemented, the current kernels are stationary and isotropic in $\mathbb{S}^2$, the surface of a sphere. Geodesic computations are done using the Vincenty formula for a sphere. The differences should be small and, honestly, and prediction error due to the spherical approximation will almost certainly be swamped by the bias and variances of small samples and uncertainty about the underlying process.
 
