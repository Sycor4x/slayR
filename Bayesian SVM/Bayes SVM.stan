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
	vector[N] idxPos;
	vector[N] idxNeg;

	//	Define some of the SVM sovler stuff once
	vector[N] alpha;
	vector[N] g;

	alpha <- rep_vector(0, N);
	g <- rep_vector(1, N);

	for(n in 1:N){
		idxPos[n] <- y[n]==1;
		idxNeg[n] <- y[n]==-1;
	}

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
	matrix[N+M,N+M] K;
	vector[N] A;
	vector[N] B;
	vector[N] yAlpha;
	int it;
	it <- 1;
	A <- rep_vector(0, N);
	B <- rep_vector(0, N);

	for(i in 1:(N+M-1)){
		for(j in (i+1):(N+M)){
			K[i,j] <- kernelFn(X_dist_sq[i,j], gamma);
			K[j,i] <- K[i,j];
		}
	}
	for(k in 1:(M+N)) K[k,k] <- 1;

	for(n in 1:N){
		if(idxNeg[n]){
			A[n] <- -C;
		}
		if(idxPos[n]){
			B[n] <- C;
		}
	}
	while(it < ITMAX & yG[i]-yG[j]<TOL){
		int i;
		int j;

		yAlpha <- alpha;
		yAlpha[idxNeg] <- -yAlpha[idxNeg];
		yG <- g;
		yG[idxNeg] <- -yG[idxNeg];

		iUp <- yAlpha < B;
		iDown <- find(yAlpha > A);

		i <- iUp[argmax(yG[iUp])];

		Ki <- K[i];
		
		idx <- iDown[yG[i] > yG[iDown]];
		Kk <- ???;
		j <- idx[argmax(pow(yG[i]-yG[idx],2)/(Ki[i]-2*Ki[idx]+Kk))]

		if(yG[i]-yG[j]<TOL)
	}
}
model{
	C ~ normal(0,1);
	gamma ~ normal(0,1);
}
generated quantities{
	
}
