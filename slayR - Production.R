##	Input:	function, bounding box for parameters
##			Control parameters: number of initialization points, number of iterations (optional: alternative criteria for convergence), number of candidate points tested at each iteration, boolean parallel flag (optional: number of parallel processors to use), 
##	Output:	the estimate of the global function minimum, a probability assessment that this is the function minimum and optionally an estimate of variance in the function v alue at the identified value. (If the function maximum is desired, just make a new objective function that multiplies the original objective function by -1.)
##	Argument box is a Dx2 matrix of lower and upper bounds on parameters. The first column is the lower bound on dimension d for 1<=d<=D, named "min". The second column is the upper bound on the same, named "max". Optionally, the rows may also be named.

library(foreach)
library(doParallel)
logit <- function(x) log(x)-log(1-x)

hart3 <- function(xx, stoch=TRUE){
  ##########################################################################
  #
  # HARTMANN 3-DIMENSIONAL FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, x3)
  #
  ##########################################################################
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)
  A <- matrix(A, 4, 3, byrow=TRUE)
  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)
  P <- matrix(P, 4, 3, byrow=TRUE)
	
  xxmat <- matrix(rep(xx,times=4), 4, 3, byrow=TRUE)
  inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
  outer <- sum(alpha * exp(-inner))
	
  if(stoch){
      y <- -outer+rnorm(10)
      return(data.frame(y_mean=mean(y), y_var=var(y)))
    }else{
      y <- -outer
      return(data.frame(y_mean=y))
    }
}

hart6 <- function(xx, stoch=TRUE){
  ##########################################################################
  #
  # HARTMANN 6-DIMENSIONAL FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, x3, x4, x5, x6)
  #
  ##########################################################################
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)
  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)
  
  xxmat <- matrix(rep(xx,times=4), 4, 6, byrow=TRUE)
  inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
  outer <- sum(alpha * exp(-inner))
  
  if(stoch){
    y <- -(2.58 + outer) / 1.94 + rnorm(10)
    return(data.frame(y_mean=mean(y), y_var=var(y)))
  }else{
    y <- -(2.58 + outer) / 1.94
    return(data.frame(y_mean=y))
  }
}

braninmodif <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi), stoch=TRUE){
  ##########################################################################
  #
  # BRANIN FUNCTION, MODIFIED
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2)
  # a = constant (optional), with default value 1
  # b = constant (optional), with default value 5.1/(4*pi^2)
  # c = constant (optional), with default value 5/pi
  # r = constant (optional), with default value 6
  # s = constant (optional), with default value 10
  # t = constant (optional), with default value 1/(8*pi)
  #
  #	usually evaluated [-5,10] x [0,15]
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  if(stoch){
    y <- term1 + term2 + s + 5*x1 + 16.64402+ rnorm(10, sd=0.5)
    return(data.frame(y_mean=mean(y), y_var=var(y)))
  }else{
    y <- term1 + term2 + s + 5*x1 + 16.64402
    return(y)
  }
}

forretal08 <- function(x, stoch=TRUE){
  ##########################################################################
  #
  # FORRESTER ET AL. (2008) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##
  ##
  ## Global min has value f(x)-6.02074 and is located at x=0.7572491
  ##
  ##
  ##########################################################################
  
  fact1 <- (6*x - 2)^2
  fact2 <- sin(12*x - 4)
 
 if(stoch){
    y <- fact1 * fact2 + rnorm(10)
    return(data.frame(y_mean=mean(y), y_var=var(y)))
  }else{
    y <- fact1 * fact2
    return(data.frame(y_mean=mean(y)))
  }
}

grlee12 <- function(x, stoch=TRUE){
  ##########################################################################
  #
  # GRAMACY & LEE (2012) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4
  
  if(stoch){
    y <- term1 + term2 + rnorm(10, sd=sqrt(0.333))
    return(data.frame(y_mean=mean(y), y_var=var(y)))
  }else{ 
    y <- term1 + term2
    return(data.frame(y_mean=y))
  }
}

goldpr <- function(xx, stoch=TRUE){
  ##########################################################################
  #
  # GOLDSTEIN-PRICE FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1a <- (x1 + x2 + 1)^2
  fact1b <- 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2
  fact1 <- 1 + fact1a*fact1b
	
  fact2a <- (2*x1 - 3*x2)^2
  fact2b <- 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2
  fact2 <- 30 + fact2a*fact2b

	if(stoch){
    y <- log(fact1*fact2)-log(3)+rnorm(10)
    return(data.frame(y_mean=mean(y), y_var=var(y)))
  }else{
    y <- log(fact1*fact2)-log(3)
    return(data.frame(y_mean=y))      
  }
}

sampleLHS <- function(size, hyperBox, D){
  library(lhs)
  boxMin <- hyperBox$min
  boxMax <- hyperBox$max
  LHS <- randomLHS(size, D)
  for(i in 1:D){
    LHS[,i] <- boxMin[i]+LHS[,i]*(boxMax[i]-boxMin[i])
  }
  LHS <- matrix(LHS, ncol=D)
  colnames(LHS) <- rownames(hyperBox)
  return(LHS)
}

applyKernel <- function(newX, FUN, ...){
  foreach(i=1:nrow(newX))%dopar%FUN(newX[i,], ...)
}

anisotropicRBF <- "
  data{
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> D;
    matrix[N,D] x;
    matrix[M,D] x_tilde;
    vector[N] y;
    vector<lower=0>[N] y_var;
  }
  transformed data{
    real y_min;
    y_min <- min(y);
  }
  parameters{
    real<lower=0, upper=100> eta_sq;
    vector<lower=0, upper=100>[D] inv_rho_sq;
    real<lower=0> sigma_sq;
    real mu;
  }
  transformed parameters{
    matrix[N,N] L;
    matrix[D,D] rho_sq_mat;

    {
      matrix[N,N] Sigma;            
      vector[D] rho_sq;
      real dist_sq;
      for(d in 1:D) rho_sq[d] <- inv(inv_rho_sq[d]);
      rho_sq_mat <- diag_matrix(rho_sq);
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          dist_sq <- quad_form(rho_sq_mat, (x[i]-x[j])' );
          Sigma[i,j] <- eta_sq * exp( -dist_sq );
          Sigma[j,i] <- Sigma[i,j];
        }
      }
      for(k in 1:N) Sigma[k,k] <- eta_sq+y_var[k]+sigma_sq;
      L <- cholesky_decompose(Sigma);
    }
  }   
  model{
    vector[N] mu_vec;
    mu_vec <- rep_vector(mu, N);
    eta_sq ~ student_t(3,0,1);
    inv_rho_sq ~ student_t(3,0,1);
    y ~ multi_normal_cholesky(mu_vec, L);
    sigma_sq ~ normal(0,1);
  }
  generated quantities{
    vector[M] Sigma_diag;
    vector[M] y_hat;    
    vector[M] delta;
    { 
      matrix[N,M] X;
      matrix[M,N] A;
      matrix[N,M] B;
      real dist_sq;
      
      for(i in 1:N){
        for(j in 1:M){
          dist_sq <- quad_form(rho_sq_mat, (x[i]-x_tilde[j])' );
          B[i,j] <- eta_sq*exp(-dist_sq);
        }
      }
      X <- mdivide_left_tri_low(L, B);
      A <- mdivide_right_tri_low(X', L);
      Sigma_diag <- eta_sq - diagonal(crossprod(X));
      y_hat <- mu + A*(y-mu);
      delta <- y_min-y_hat;
    }
  }
"

anisotropicMatern <- "
  data{
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> D;
    matrix[N,D] x;
    matrix[M,D] x_tilde;
    vector[N] y;
    vector<lower=0>[N] y_var;
    int<lower=0> p;
  }
  transformed data{
    real y_min;
    real q;
    real<lower=0.5> nu;
    y_min <- min(y);
    q <- lgamma(p+1)-lgamma(2*p+1);
    nu <- p+0.5;
  }
  parameters{
    real<lower=0, upper=100> eta_sq;
    vector<lower=0, upper=100>[D] inv_rho_sq;
    real<lower=0> sigma_sq;
    real mu;
  }
  transformed parameters{
    matrix[N,N] L;
    matrix[D,D] Lambda;
    {
      matrix[N,N] Sigma;            
      vector[D] rho_sq;
      real dist;
      vector[p+1] fact2;

      for(d in 1:D) rho_sq[d] <- inv(inv_rho_sq[d]);
      Lambda <- 2*nu*diag_matrix(rho_sq);
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          dist <- sqrt(quad_form(Lambda, (x[i]-x[j])' ) );
          for(k in 0:p) fact2[k+1] <- exp( lgamma(k+p+1) - lgamma(k+1) - lgamma(p-k+1) ) * (2*dist)^(p-k);
          Sigma[i,j] <- eta_sq * exp(q - dist) * sum(fact2);
          Sigma[j,i] <- Sigma[i,j];
        }
      }
      for(k in 1:N) Sigma[k,k] <- eta_sq+y_var[k]+sigma_sq;
      L <- cholesky_decompose(Sigma);
    }
  }   
  model{
    vector[N] mu_vec;
    mu_vec <- rep_vector(mu, N);
    sigma_sq ~ normal(0,1);
    eta_sq ~ student_t(3,0,1);
    inv_rho_sq ~ student_t(3,0,1);
    y ~ multi_normal_cholesky(mu_vec, L);
  }
  generated quantities{
    vector[M] Sigma_diag;
    vector[M] y_hat;    
    vector[M] delta;
    { 
      matrix[N,M] X;
      matrix[M,N] A;
      matrix[N,M] B;
      real dist;
      vector[p+1] fact2;
      
      for(i in 1:N){
        for(j in 1:M){
          dist <- sqrt(quad_form(Lambda, (x[i]-x_tilde[j])' ));
          for(k in 0:p) fact2[k+1] <- exp(lgamma(k+p+1) - lgamma(k+1) - lgamma(p-k+1))*(2*dist)^(p-k);
          B[i,j] <- eta_sq * exp(q-dist) * sum(fact2);
        }
      }
      X <- mdivide_left_tri_low(L, B);
      A <- mdivide_right_tri_low(X', L);
      Sigma_diag <- eta_sq - diagonal(crossprod(X));
      y_hat <- mu + A*(y-mu);
      delta <- y_min-y_hat;
    }
  }
"

library(rstan)
anisotropicRBF_stan_compiled <- stan(model_code=anisotropicRBF_stan, iter=10)

nextPointLCB <- function(MCsimulation, newStarts=startNewX, iter, D, delta, earlyStop_in, prev_metric, epsilon, ...){
  y_hat <- MCsimulation$y_hat
  Sigma_diag <- MCsimulation$Sigma_diag
  s <- ifelse(Sigma_diag>0, sqrt(Sigma_diag), 0)
  beta_iter <- 2*(log(D*(pi*iter)^2)-log(6*delta))
  LCB <- sapply(1:ncol(y_hat), function(i) y_hat[,i]-sqrt(beta_iter)*s)
  LCB <- colMeans(LCB)
  minimum <- which(LCB==min(LCB))
  if(length(minimum)>1){
    minimum <- minimum[1]
    print(paste0("Informational message: LCB acquisition function reports ", length(minimum), " minima. Selecting one minimum at random and proceeding.")) #Note that the points are already random samples, so selecting the first is the same as selecting randomly.
  }
  new_x <- newStarts[minimum,]
  new_metric <- LCB[minimum]
  if(earlyStop_in & any(!is.na(prev_metric))){
    terminate <- new_metric > min(prev_metric, na.rm=TRUE)*epsilon
    ## We want to terminate when new_metric < min(prev_metric)*epsilon. So terminate is set to FALSE otherwise, implying that we  need to check new_metric > min(prev_metric, na.rm=TRUE)*epsilon)
  }else{
    terminate <- FALSE
  }
  return(list(new_x=new_x, new_metric=new_metric, terminate=terminate))
}

nextPointEI <- function(MCsimulation, newStarts=startNewX, iter, D, delta, earlyStop_in, prev_metric, epsilon, ...){
  delta <- MCsimulation$delta
  y_hat <- MCsimulation$y_hat
  Sigma_diag <- MCsimulation$Sigma_diag
  s <- ifelse(Sigma_diag>0, sqrt(Sigma_diag), 0)
  EI <- sapply(1:ncol(y_hat), function(i) delta[,i]*pnorm(delta[,i], sd=s[,i])+s[,i]*dnorm(delta[,i], sd=s[,i]))
  EI <- colMeans(EI)
  maximum <- which(EI==max(EI))
  if(length(maximum)>1){
    maximum <- maximum[1]
    print(paste0("Informational message: EI acquisition function reports ", length(maximum), " maxima. Selecting one maximum at random and proceeding.")) #Note that the points are already random samples, so selecting the first is the same as selecting randomly.
  }
  new_x <- newStarts[maximum,]
  new_metric <- EI[maximum]
  if(new_metric<0) print(paste("Non-fatal computational hiccup (really!). Largest expected improvement is", new_metric, "which is negative, but all EI values (not just the maximum) are defined to always be non-negative. Plausible sources of error include (1) accumulated numerical error and (2) a bug in the GP program. Computation will proceed, ignoring the negative EI result."))
  if(earlyStop_in & any(!is.na(prev_metric))){
    terminate <- new_metric < epsilon & new_metric > 0
  }else{
    terminate <- FALSE
  }
  return(list(new_x=new_x, new_metric=new_metric, terminate=terminate))
}

#test1 <- matrix(rnorm(30), nrow=10, ncol=3)
#test2 <- rnorm(7)
#anisotropicRBF(mat1=test1, mat2=test2, par_list=list(eta_sq=1.5, rho_sq_mat=diag(rnorm(3))))

varPrep <- function(fnReturn, flag, epsilon, y_mean_in=y_mean, y_var_in=y_var){
	if(flag){
		y_mean <- c(y_mean_in, fnReturn$y_mean)
		y_var <- c(y_var_in, fnReturn$y_var)
		y_var[y_var<epsilon] <- epsilon
	}else{
		y_mean <- c(y_mean_in, fnReturn$y_mean)
		y_var	 <- rep(epsilon, length(y_mean))
	}
	return(list(y_mean=y_mean, y_var=y_var))
}

##	function f is the function to be minimized. If you want f to be maximized, just flip the sign of the output. Function f also optionally takes additional arguments. These are supplied to EGO via fOptions.

slayR <- function(f, box, initial, iter, nStan=200, nCandidate=200, isStochastic=TRUE, epsDiag=1e-6, kernel="Gaussian", verbose=FALSE, plotProgress=TRUE, earlyStop=FALSE, earlyStopTol, acquisitionMethod="LCB", LCB_delta, Matern_p, ...){
  ## Intializing options
  stanProgram <- switch(kernel, 
    Gaussian = anisotropicRBF, 
    Matern = anisotropicMatern)
  acquisitionFn <- switch(acquisitionMethod, 
    LCB = nextPointLCB, 
    EI = nextPointEI)
  dataPackageExpression <- switch(kernel, 
    Gaussian=expression(list(N=initial+i-1, M=nCandidate, D=D, x=x, y=y_mean, y_var=y_var, x_tilde=startNewX)), 
    Matern=expression(list(N=initial+i-1, M=nCandidate, D=D, x=x, y=y_mean, y_var=y_var, x_tilde=startNewX, p=Matern_p)))

  addlArgs <- as.list(substitute(list(...)))[-1L]

  ## The Novetta color pallette
  Ndarkblue <- rgb(0, 100, 182, maxColorValue=255)
  Ngray <- rgb(147, 149, 152, maxColorValue=255)
  Nteal  <- rgb(0, 171, 194, maxColorValue=255)
  Norange <- rgb(247, 155, 46, maxColorValue=255)

  ## More initializing
  D <- nrow(box)
  if(missing(LCB_delta) & acquisitionMethod=="LCB"){ 
    LCB_delta <- 0.1
    print(paste0("Informational message: acquisitionMethod='LCB' but parameter LCB_delta is missing. Setting it to default value LCB_delta=", LCB_delta, "."), quote=FALSE)
  }
  if(acquisitionMethod=="LCB" & (LCB_delta < 0 | LCB_delta > 1) | length(LCB_delta)!=1) {
    stop(paste("User error: Either length(LCB_Delta)!=1, or LCB_delta outside of bounds: 0 < LCB_delta < 1. User-provided value of LCB_delta:", LCB_delta))
  }
  if(earlyStop & missing(earlyStopTol) & acquisitionMethod%in%(c("LCB", "EI"))){
    earlyStopTol <- switch(acquisitionMethod, LCB = 1, EI = 5e-5)
    print(paste0("Informational message: earlyStopTol missing and earlyStop=TRUE. Setting default tolerance for acquisitionMethod=", acquisitionMethod, ": earlyStopTol=", earlyStopTol), quote=FALSE)
  }
  if(kernel=="Matern"){
    if(missing(Matern_p)){
      Matern_p <- 2L
      print(paste0("Informational message: kerenel='Matern' but Matern_p not provided. Setting Matern_p to default value of ", Matern_p, "."), quote=FALSE)
    }
    if(Matern_p<0 | !is.integer(Matern_p)) stop("User errror: kernel='Matern' set and Matern_p is provided by user, but its value is invalid. When kernel='Matern' is set, argument Matern_p must be a non-negative integer. A future release may admit non-integer values Matern_p. Default value Matern_p=2L.")
  }
	if(!is.data.frame(box)|ncol(box)!=2|nrow(box)<1|any(names(box)!=c("min", "max"))) stop("User error: box must be a Dx2 data.frame containing the boundary limits for each dimension under optimization. Also, names(box)=c('min', 'max').", quote=FALSE)
  if(earlyStop & !(acquisitionMethod%in%c("LCB","EI"))){
    stop("User error: earlyStop may only be set to TRUE if used with acquisitionMethod='LCB' or acquisitionMethod='EI'. If acquisitionMethod is set to something else, you must set earlyStop=FALSE. This may be revised in a future release.")
    }
	if(missing(initial)) initial <- min(20*D, 100)
	if(missing(iter)) iter <- min(250, 20*D)

	# This just evaluates the sample points at the derived values
	if(is.list(initial)){
		print("Checking if provided intput can be used as starting values...", quote=FALSE)
		print("Informational message: if initial is a list, must supply named components x and fOut. Further, fOut must have named component y_mean. If the isStochastic flag is TRUE, fOut must further have named component y_var. If initial is not a list, then set to NULL (default) or a positive integer to initialize randomly using Latin hyper-cube sampling.", quote=FALSE)
		if(any(!(names(initial)%in%c("x", "fOut")))) stop("User error: if initial is a list, must supply named components x and fOut.", quote=FALSE)
		x <- initial$x
		fOut <- initial$fOut
		initial <- nrow(x)
	}else{
		print(paste0("Informational message: default initialization underway using random sampling of ", initial, " initial design points using Latin hypercubes."), quote=FALSE)
		x <- matrix(sampleLHS(initial, hyperBox=box, D), nrow=initial, ncol=D, dimnames=list(NULL, rownames(box)))
		fOut <- apply(x, MARGIN=1, FUN=f, ...)
		fOut <- do.call(rbind, fOut)
	}

	metric <- rep(NA, initial)
	postProc <- varPrep(fOut, flag=isStochastic, epsilon=epsDiag, y_mean_in=NULL, y_var_in=NULL)
	y_mean <- postProc$y_mean
	y_var <- postProc$y_var
	colnames(x) <- rownames(box)
  if(verbose){
    print("Initial design points:", quote=FALSE)
    print(cbind(x,y_mean,y_var))
  }

  startNewX <- sampleLHS(nCandidate, hyperBox=box, D=D)
  i <-1
  stanData <- eval(dataPackageExpression)
  if(verbose) print("Informational message: Compiling stan program now...")
  GPmodel0 <- stan(model_code=stanProgram, iter=nStan, data=stanData, cores=4)
  if(verbose) print("Informational message: Compilation complete! Beginning inference.")

	for(i in 1:iter){
		print(paste0("Iteration ", i, " of at most ", iter, "."), quote=FALSE)
		#(1) collect data for stan estimation of GP model; estimate GP model 
		startNewX <- sampleLHS(nCandidate, hyperBox=box, D=D)
    stanData <- eval(dataPackageExpression)
		GPmodel <- stan(fit=GPmodel0, data=stanData, iter=nStan, cores=4)
		GPmodelPars <- extract(GPmodel, pars=c("y_hat", "delta", "Sigma_diag"), inc_warmup=FALSE)
    #(2) evaluate acquistion function to select next point
		acqOut <- acquisitionFn(MCsimulation=GPmodelPars, newStarts=startNewX, D=D, iter=i, delta=LCB_delta, earlyStop_in=earlyStop, epsilon=earlyStopTol, prev_metric=metric)
    new_metric <- acqOut$new_metric
    new_x <- acqOut$new_x
    terminate <- acqOut$terminate
		#(3) evaluate function at new sample point
		fOut <- f(new_x, ...)
		#(4) append results, return to (1)
		postProc <- varPrep(fOut, flag=isStochastic, epsilon=epsDiag, y_mean_in=y_mean, y_var_in=y_var)
		y_mean <- postProc$y_mean
		y_var <- postProc$y_var
		x <- rbind(x, new_x)
		metric <- c(metric, new_metric)
    #(5) Interface options: plotting progress, printing progress, checking stopping criteria.
    if(verbose){
      print("New point acquired:", quote=FALSE)
      print(cbind(x,y_mean,y_var,metric)[initial+i,])
    }
    if(plotProgress){
      par(mfrow=c(2,1), mar=c(3, 4, 4, 1))
      plot(1:i, y_mean[initial+(1:i)], main="Progress", ylab="y", xlim=c(1, iter), type="b", col=Ndarkblue)
      abline(h=min(y_mean), col=Norange, lwd=2)
      logFlag <- ifelse(all(metric[!is.na(metric)]>0) & acquisitionMethod=="EI", "y", "")
      plot(1:i, metric[initial+(1:i)], xlim=c(1,iter), type="b", ylab=acquisitionMethod, col=Ndarkblue, log=logFlag, xlab="Iteration")
      if(earlyStop & acquisitionMethod=="EI" & !missing(earlyStopTol)) if(earlyStopTol>0) abline(h=earlyStopTol, col=Norange)
    }
    if(earlyStop & terminate){
      print("Informational message: Early stopping criteria satisfied. Exiting optimizaiton now.", quote=FALSE)
      break
    }
    if(i >= iter) print("Maximum number of iterations met. Exiting optimization now.", quote=FALSE)	
  }
	## Optimization complete; collect results and return to user
	rownames(x) <- c(rep("Initial", initial), paste0("Iteration", 1:(nrow(x)-initial)))
	xfx <- cbind(x, y_mean, y_var, metric)
  xfx_iter <- xfx[rownames(xfx)!="Initial",]
  xfx_init <- xfx[rownames(xfx)=="Initial",]
	
	best_ndx <- which(y_mean==min(y_mean))
	value <- y_mean[best_ndx]
	best <- x[best_ndx,]
	out <- list(xfx=xfx, xfx_init=xfx_init, opt_path=xfx_iter, best=best, value=value, new_metric=new_metric, iter=i, GPfinal=GPmodel)
	print(paste0("Exiting slayR. Best function value obtained: ", value, ". Acquisition metric at final iteration: ", new_metric), quote=FALSE)
  dev.off()
	return(out)
}

# Goldstein-Price function
##  Default acquisition fn
goldprBox <- data.frame(min=c(-2,-2),max=c(2,2))
rownames(goldprBox) <- c("x1", "x2")
foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, iter=10, earlyStop=FALSE, kernel="Matern")

foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, iter=20, earlyStop=FALSE, kernel="Gaussian")

foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, iter=20, earlyStop=FALSE)

foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, iter=10, earlyStop=FALSE, kernel="Matern")

repackage <- list(x=foo$xfx[,(1:2)], fOut=list(y_mean=foo$xfx[,3], y_var=foo$xfx[,4]))
foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, initial=repackage, iter=20)

repackage <- list(x=foo$xfx[,(1:2)], fOut=list(y_mean=foo$xfx[,3]))
foo <- slayR(goldpr, box=goldprBox, isStochastic=FALSE, initial=repackage, iter=10, stanProgram=anisotropicRBF_stan_sigma_est, epsDiag=0, EGOtol=0)

# Goldstein-Price function
##  Expected improvement acquisition fn
goldprBox <- data.frame(min=c(-2,-2),max=c(2,2))
rownames(goldprBox) <- c("x1", "x2")
foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, iter=10, earlyStop=FALSE, acquisitionMethod="EI")

repackage <- list(x=foo$xfx[,(1:2)], fOut=list(y_mean=foo$xfx[,3], y_var=foo$xfx[,4]))
foo <- slayR(goldpr, box=goldprBox, isStochastic=TRUE, initial=repackage, iter=10, acquisitionMethod="EI")

library(fields)
X1 <- seq(-2,2,len=99)
X2 <- seq(-2,2,len=101)
grid <- expand.grid(X1=X1, X2=X2)
gridOut <- apply(grid,1, goldpr, stoch=FALSE)
gridOut <- unlist(gridOut)
gridOut <- data.frame(grid, Y=gridOut)
matFn <- as.matrix(unstack(gridOut, Y~X2))
image.plot(matFn, x=X1, y=X2, xlab=expression(X[1]), ylab=expression(X[2]), col=tim.colors(250), zlim=c(0,13), main="Goldstein-Price Function, log scale")
points(foo$opt_path[,(1:2)])
lines(foo$opt_path[,(1:2)])
points(x=foo$best[1], y=foo$best[2], lwd=3, col="black")
points(foo$xfx_init[,(1:2)], col="grey", lwd=2)
points(x=0, y=-1, col="white", lwd=2)

#######################
##	Modified Branin-Hoo function
##	

braninhooBox <- data.frame(min=c(-5,0), max=c(10, 15))
rownames(braninhooBox) <- c("x1", "x2")
foo <- slayR(braninmodif, box=braninhooBox, isStochastic=TRUE, iter=10)

library(fields)
X1 <- seq(-5,10,len=101)
X2 <- seq(-0,15,len=99)
grid <- expand.grid(X1=X1, X2=X2)
gridOut <- apply(grid,1, braninmodif)
gridOut <- lapply(gridOut, function(x) x$y_mean)
gridOut <- data.frame(grid, Y=do.call(rbind, gridOut))
matFn <- as.matrix(unstack(gridOut, Y~X2))
image.plot(matFn, x=X1, y=X2, xlab=expression(X[1]), ylab=expression(X[2]), col=tim.colors(250))
points(foo$opt_path[,(1:2)])
lines(foo$opt_path[,(1:2)])
points(x=foo$best[1], y=foo$best[2], lwd=3, col="black")
points(x=-3.689242, y=13.630397, col="white", lwd=3)

##	Deterministic grlee12 function

foo <- slayR(grlee12, stoch=FALSE, box=data.frame(min=0.5, max=2.5), isStochastic=FALSE, iter=10)
x <- seq(0.5, 2.5, by=0.001)
plot(x, unlist(grlee12(x, stoch=FALSE)), type="l", ylab="Gramacy & Lee 2012 function", xlab="x")
points(foo$xfx_init[,1], foo$xfx_init[,2])
text(foo$opt_path[,1], foo$opt_path[,2], 1:foo$iter, col="black", lwd=2)
points(foo$best, foo$value, col="red", lwd=2)

##	Stoch forretal08 function
##	
forretal08box <- data.frame(min=c(0), max=c(1))
rownames(forretal08box) <- "x"
foo <- slayR(forretal08, box=forretal08box, isStochastic=FALSE, iter=10, stanProgram=anisotropicRBF_stan_sigma_est)

#################
## HARTMANN 3 function
hart3box <- data.frame(min=c(0,0,0), max=c(1,1,1))
rownames(hart3box) <- c("x1", "x2", "x3")
foo <- slayR(hart3, box=hart3box, earlyStop=TRUE)

#######################
##	HARTMANN 6 function
hart6box <- data.frame(min=c(0,0,0,0,0,0), max=c(1,1,1,1,1,1))
rownames(hart6box) <- c("x1", "x2", "x3", "x4", "x5", "x6")
foo <- slayR(hart6, box=hart6box, earlyStop=TRUE)
save(foo, file="C:/Users/David Elkind/Desktop/Surrogate Models Code/Benchmarks/hartman6.RData")


