functions{
  real geodesic(row_vector x1, row_vector x2){
    vector[2] phi; // latitudes
    vector[2] lambda; // longitudes
    real delta_lambda;
    real term1;
    real term2;
    real denominator;
    real r;
    real delta_sigma;
    real dist_rho_quotientance_out;

    r <- 6378.1370; //equatorial radius of the Earth in km (under the spherical Earth model) largest error will be <1%

    phi[1] <- x1[1];
    phi[2] <- x2[1];

    lambda[1] <- x1[2];
    lambda[2] <- x2[2];

    delta_lambda <- fabs(lambda[1]-lambda[2]);

    term1 <- pow(cos(phi[2])*sin(delta_lambda),2);
    term2 <- pow(cos(phi[1])*sin(phi[2])-sin(phi[1])*cos(phi[2])*cos(delta_lambda),2);
    denominator <- sin(phi[1])*sin(phi[2])+cos(phi[1])*cos(phi[2])*cos(delta_lambda);
    delta_sigma <- atan2(sqrt(term1+term2), denominator);
    dist_rho_quotientance_out <- r*delta_sigma;
    return dist_rho_quotientance_out;
  }
}
data{
  int<lower=1> N;              //number of units in training data
  int<lower=1> M;              //number of units in test data
  matrix[N,2] x_train;         //lat-long (in that order) locations of training data X
  matrix[M,2] x_test;          //lat-long (in that order) locations of test data X~
  int<lower=0> y_train[N];     //training response, input as a vector of N binary entries
  int<lower=0> y_test[M];      //testing response, inputs as a vector of M binary entries
}
transformed data{
  matrix[N+M,N+M] x_dist;
  row_vector[3] coef;

  coef[1] <- 1;
  coef[2] <- sqrt(5);
  coef[3] <- 5.0/3.0;

  { 
    matrix[N+M,2] x_all;
    x_all <- append_row(x_train, x_test);
    
    for(i in 1:(N+M-1)){
      for(j in (i+1):(N+M)){
        x_dist[i,j] <- geodesic(x_all[i], x_all[j]);
        x_dist[j,i] <- x_dist[i,j];
      }
    }
    for(k in 1:(N+M)) x_dist[k,k] <- 0;
  } 
}
parameters{
  real<lower=0, upper=100> inv_rho;
  real<lower=0> sigma_sq;
  real mu;
  vector[N] z_train;
  vector[M] z_test;
}
transformed parameters{
  matrix[N+M,N+M] L;
  matrix[N+M,N+M] Sigma;            

  {
    real dist_rho_quotient;
    vector[3] fact2;

    for(i in 1:(N+M-1)){
      for(j in (i+1):(N+M)){
        dist_rho_quotient <- inv_rho*x_dist[i,j];
        for(k in 0:2) fact2[k+1] <- pow(dist_rho_quotient, k);
        Sigma[i,j] <- exp(-sqrt(5)*dist_rho_quotient)*(coef*fact2);
        Sigma[j,i] <- Sigma[i,j];
      }
    }
    for(k in 1:(N+M)) Sigma[k,k] <- 1+sigma_sq;
    L <- cholesky_decompose(Sigma);
  }
}   
model{
  vector[N+M] mu_vec;
  vector[N+M] z_all;
  z_all <- append_row(z_train, z_test);
  mu_vec <- rep_vector(mu, N+M);

  for (n in 1:N) y_train[n] ~ poisson_log(z_train[n]);
  z_all ~ multi_normal_cholesky(mu_vec, L);

  sigma_sq ~ chi_square(1);
  inv_rho ~ student_t(3,0,1);
}
generated quantities{
  vector[M] y_hat;
  real rmse;

  {
  vector[M] new_rmse;
    for(m in 1:M){
      y_hat[m] <- exp(z_test[m]);
      new_rmse[m] <- pow(y_test[m]-y_hat[m],2);
    }
    rmse <- mean(new_rmse);
  }
}