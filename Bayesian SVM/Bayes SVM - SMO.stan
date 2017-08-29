functions {
	real kernelFn(real gamma, real dist_sq) {
		return exp(-gamma * dist_sq);
	}
}
data {
	int<lower = 1> nTrain;
	int<lower = 1> nTest;
	int<lower = 1> nFeature;
	int<lower = 1> ITMAX;
	real<lower = 0> TOL;

	matrix[nTrain, nFeature] xTrain;
	matrix[nTest, nFeature] xTest;

	vector<lower = -1, upper = 1>[nTrain] yTrain;
	int<lower = -1, upper = 1> yTest[nTest];

//	real<lower = pow(2, -5), upper = pow(2, 12)> C;
//	real<lower = pow(2, -15), upper = pow(2, 3)> gamma;
}

transformed data {
	matrix[nTrain, nTrain] xTrainDistSq;
	matrix[nTrain, nTest] xTrainTestDistSq;
	int yTestBin[nTest];

	for(k in 1:nTest) if(yTest[k] > 0) yTestBin[k] <- 1;
	for(i in 1:nTrain) {
		xTrainDistSq[i, i] <- 0;
		for(j in 1:i - 1) {
			xTrainDistSq[i, j] <- squared_distance(xTrain[i], xTrain[j]);
			xTrainDistSq[j, i] <- xTrainDistSq[i, j];
		}
		for(j in 1:nTest) {
			xTrainTestDistSq[i, j] <- squared_distance(xTrain[i], xTest[j]);
		}
	}
}
parameters {
	vector<lower = pow(2, -5), upper = pow(2, 15)>[1] C;
	vector<lower = pow(2, -15), upper = pow(2, 3)>[1] gamma;
	vector[2] beta;
	cholesky_factor_corr[4] Omega;
	vector<lower=0>[4] tau;
	vector[4] mu;
}
transformed parameters {
	vector[nTrain] alpha;
	vector[nTest] dv;
	real rho;
	matrix[nTest,2] Z;

	{	
		matrix[2,2] Sigma;
		vector[nTrain] g;
		vector[nTrain] A;
		vector[nTrain] B;
		matrix[nTrain, nTrain] K;
	
		for(i in 1:nTrain) {
			// Initialize A, B and g
			alpha[i] <- 0;
			if(yTrain[i] == 1) {
				A[i] <- 0;
				B[i] <- C[1];
				g[i] <- 1;
			} else {
				A[i] <- -C[1];
				B[i] <- 0;
				g[i] <- -1;
			}
	
			// Precompute the entire kernel matrix
			K[i, i] <- 1;
			for(j in 1:i - 1) {
				K[i, j] <- kernelFn(gamma[1], xTrainDistSq[i, j]);
				K[j, i] <- K[i, j];
			}
		}
	
		{
			int it;
			int i;
			int j;
			int nFree;
			real minG;
			real maxG;
			real l;
			real dualGap;
			vector[3] lVector;
	
			// SMO algorithm starts here
			it <- 1;
			dualGap <- positive_infinity();
			while(it < ITMAX && dualGap > TOL) {
		   		minG <- positive_infinity();
		   		maxG <- negative_infinity();
		   		for(k in 1:nTrain) {
		   			if(alpha[k] < B[k] && g[k] > maxG) {
						i <- k;
						maxG <- g[k];
		   			}
		   			if(A[k] < alpha[k] && g[k] < minG) {
		   				j <- k;
		   				minG <- g[k];
		   			}
		   		}
	
				dualGap <- g[i] - g[j];
	
				lVector[1] <- B[i] - alpha[i];
				lVector[2] <- alpha[j] - A[j];
				lVector[3] <- dualGap / (K[i, i] + K[j, j] - 2 * K[i, j]);
				l <- min(lVector);
	
				g <- g + l * to_vector(K[j] - K[i]);
				alpha[i] <- alpha[i] + l;
				alpha[j] <- alpha[j] - l;
	
				it <- it + 1;
			}
			// Compute intercept, rho
			rho <- 0;
			nFree <- 0;
			for(k in 1:nTrain) {
				if(A[k] < alpha[k] && alpha[k] < B[k]) {
					nFree <- nFree + 1;
					rho <- rho + g[k];
				}
			}
			rho <- rho / nFree;
	
			// Compute predictions on test data
			for(m in 1:nTest) {
				dv[m] <- rho;
				for(n in 1:nTrain) {
					dv[m] <- dv[m] + alpha[n] * kernelFn(gamma[1], xTrainTestDistSq[n, m]);
				}
			}
		}
	}
	Z <- append_col(rep_vector(1, nTest), dv);
}
model {
	vector[4] hyper_pars;
	hyper_pars <- append_row(append_row(beta,C), gamma);
	yTestBin ~ bernoulli_logit(Z*beta);
	hyper_pars ~ multi_normal_cholesky(mu, diag_matrix(tau)*Omega);
	tau[1:2] ~ cauchy(0,2.5);
	Omega ~ lkj_corr_cholesky(2.0);
}
generated quantities {
}
