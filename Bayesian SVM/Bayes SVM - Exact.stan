functions{
	real kernelFn(real dist_sq, real par){
		real out;
		out <- exp(-par*dist_sq);
		return out;
	}
}
data{
	int<lower=1> N;
	int<lower=1> M;
	int<lower=1> D;
	int<lower=10> ITMAX;

	matrix[N,D] X;
	matrix[M,D] X_tilde;

	int<lower=-1, upper=1> y[N];
	int<lower=-1, upper=1> y_tilde[M];
}
transformed data{
	// Transformations of the data
	matrix[N+M,N+M] X_dist_sq;

	//	Define some of the SVM sovler stuff once
	{
		matrix[N+M,D] X_all;
		X_all <- append_row(X, X_tilde);
		for(i in 1:(N+M-1)){
			for(j in (i+1):(N+M)){
				X_dist_sq[i,j] <- squared_distance(X_all[i], X_all[j]);
				X_dist_sq[j,i] <- X_dist_sq[i,j];
			}
		}
	}
	for(k in 1:(N+M)) X_dist_sq[k,k] <- 0;
}
parameters{
	real<lower=0> C;
	real<lower=0> gamma;
}
transformed parameters{
	vector[N] alpha;
	matrix[N,N] P;
	matrix[N+M,N+M] K;
	vector[N] idx1;
	vector[N] idx2;
	vector[N] e;
	matrix[N,N] I;

	alpha <- rep_vector(0, N);
	e <- rep_vector(1, N);
	I <- diag_matrix(1);
	idx1 <- alpha < 1e-8;
	idx2 <- alpha > C-1e-8;

	for(i in 1:(N+M-1)){
		for(j in (i+1):(N+M)){
			K[i,j] <- kernelFn(X_dist_sq[i,j], gamma);
			K[j,i] <- K[i,j];
		}
	}
	for(k in 1:(N+M)) K[k,k] <- 1;


	
	{ //Begin SVM solver
	}
	{ //Begin SVM predictions

	}
}
model{
	C ~ normal(0,1);
	gamma ~ normal(0,1);
}
generated quantities{
	
}
