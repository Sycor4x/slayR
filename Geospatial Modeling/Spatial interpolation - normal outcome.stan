functions{   
  real Matern(real distance, real inv_rho){
    real out;
    row_vector[3] coef;
    vector[3] fact2;
    coef[1] <- 1;
    coef[2] <- sqrt(5);
    coef[3] <- 5.0/3.0;

    for(i in 0:2) fact2[i+1] <- pow(inv_rho*distance, i);
    out <- exp(-sqrt(5)*inv_rho*distance)*(coef*fact2);
    return out;
  }
  real geodesic(row_vector x1, row_vector x2){
    vector[2] phi;
    vector[2] lambda; 
    real delta_lambda;
    real term1;
    real term2;     
    real denominator;
    real r;
    real delta_sigma;
    real dist_out;

    r <- 6378.1370; //equatorial radius of the Earth in km (under the spherical Earth model) largest error will be <1%

    phi[1] <- x1[1]*pi()/180;
    phi[2] <- x2[1]*pi()/180;

    lambda[1] <- x1[2]*pi()/180;
    lambda[2] <- x2[2]*pi()/180;

    delta_lambda <- fabs(lambda[1]-lambda[2]);

    term1 <- pow(cos(phi[2])*sin(delta_lambda),2);
    term2 <- pow(cos(phi[1])*sin(phi[2])-sin(phi[1])*cos(phi[2])*cos(delta_lambda),2);
    denominator <- sin(phi[1])*sin(phi[2])+cos(phi[1])*cos(phi[2])*cos(delta_lambda);
    delta_sigma <- atan2(sqrt(term1+term2), denominator);
    dist_out <- r*delta_sigma;
    return dist_out;
  }
}
data{
  int<lower=1> N;              //number of units in training data
  int<lower=1> M;              //number of units in test data
  matrix[N,2] x_train;         //lat-long (in that order) locations of training data X
  matrix[M,2] x_test;          //lat-long (in that order) locations of test data X~
  vector[N] y_train;           //training response, input as a vector of N entries
  vector[M] y_test;            //testing response, inputs as a vector of M entries   
}
transformed data{
  matrix[N,N] x_dist;
  matrix[N,M] x_tilde_dist;
  
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      x_dist[i,j] <- geodesic(x_train[i], x_train[j]);
      x_dist[j,i] <- x_dist[i,j];
    }
  }
  for(k in 1:N) x_dist[k,k] <- 0;
  for(i in 1:N){
    for(j in 1:M){
      x_tilde_dist[i,j] <- geodesic(x_train[i], x_test[j]);
    }
  }
}
parameters{
  real<lower=0, upper=100> inv_rho;
  real<lower=0> sigma_sq;
  real mu;
}
transformed parameters{
  matrix[N,N] L;

  {
    matrix[N,N] Sigma;
    Sigma <- diag_matrix(rep_vector(1+sigma_sq, N));
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        Sigma[i,j] <- Matern(x_dist[i,j], inv_rho);
        Sigma[j,i] <- Sigma[i,j];
      }
    }
    L <- cholesky_decompose(Sigma);
  }
}
model{
  vector[N] mu_vec;
  mu_vec <- rep_vector(mu, N);
  y_train ~ multi_normal_cholesky(mu_vec, L);
  sigma_sq ~ chi_square(1);
  inv_rho ~ student_t(3,0,1);
}
generated quantities{
  vector[M] Sigma_diag;
  vector[M] y_hat;
  real rmse;
  { 
    matrix[N,M] X;
    matrix[M,N] A;
    matrix[N,M] B;
    
    for(i in 1:N){
      for(j in 1:M){
        B[i,j] <- Matern(x_tilde_dist[i,j], inv_rho);
      }
    }
    X <- mdivide_left_tri_low(L, B);
    A <- mdivide_right_tri_low(X', L);
    Sigma_diag <- 1 - diagonal(crossprod(X));
    y_hat <- mu + A*(y_train-mu);
  }
  {
    vector[M] junk;
    for(k in 1:M) junk[k] <- pow(y_hat[k]-y_test[k],2);
    rmse <- sqrt(mean(junk));
  }
}