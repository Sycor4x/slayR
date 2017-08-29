#Bayesian SVM
library(rstan)
library(data.table)

X <- data.table(iris)
X[,Species:=ifelse(X[,Species]=="versicolor", 1, -1)]
Y <- X[,Species]
X <- as.matrix(X[,Species:=NULL])
X <- as.matrix(scale(X))

test <- 1:nrow(X)
# test <- sample(nrow(X), size=nrow(X)/5, replace=FALSE)
# train <- (1:nrow(X))[-test]

train <- test

C <- 100
gamma <- 1.208

iris_test <- list(nTrain=length(train), nTest=length(test), nFeature=ncol(X), ITMAX=1e5, TOL=1e-8, xTrain=X[train,], yTrain=Y[train], xTest=X[test,], yTest=Y[test], C=C, gamma=gamma)

stanOut <- stan(file="C:/Users/David Elkind/Desktop/RConvFn/Bayesian SVM/Bayes SVM - SMO.stan", data=iris_test, chains=1, iter=10)
alpha <- extract(stanOut, "alpha")[[1]][1,]
extract(stanOut, "obj_fn")[[1]]

K <- extract(stanOut, "K")[[1]][1,,]

stanObj <- 1/2*t(alpha)%*%K%*%alpha-Y[train]%*%alpha

libsvmObj <- 1/2*t(libsvm_alpha)%*%K%*%libsvm_alpha-Y[train]%*%libsvm_alpha

(stanObj-libsvmObj)

(alpha-libsvm_alpha)%*%(alpha-libsvm_alpha)

extract(stanOut, "attained_it")[[1]]

library(e1071)
libsvm <- svm(y=as.factor(Y[train]), x=X[train,], gamma=gamma, cost=C, scale=FALSE, shrinking=FALSE, tolerance = 1e-8)

libsvm_alpha <- rep(0, length(train))
libsvm_alpha[libsvm$index] <- libsvm$coefs*Y[train][1]

length(libsvm_alpha)==length(alpha)

abs_dif <- (alpha-libsvm_alpha)
plot(abs_dif)
abline(h=0)

rmse_alpha <- sqrt((alpha-libsvm_alpha)%*%(alpha-libsvm_alpha))
rmse_alpha

t(libsvm_alpha)%*%K%*%libsvm_alpha/2
