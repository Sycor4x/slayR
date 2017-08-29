#Bayesian SVM
library(rstan)
library(data.table)

# newIris <- read.csv("C:/Users/David Elkind/Desktop/RConvFn/expanded iris data/newIris.csv", stringsAsFactors=FALSE)
# X <- data.table(newIris)

X <- data.table(iris)

X[,Species:=ifelse(X[,Species]=="versicolor", 1, -1)]
Y <- X[,Species]
X <- as.matrix(X[,Species:=NULL])
X <- as.matrix(scale(X))

test <- sample(nrow(X), size=nrow(X)/2, replace=FALSE)
train <- (1:nrow(X))[-test]

iris_test <- list(nTest=length(test), nTrain=length(train), nFeature=ncol(X), ITMAX=1000, TOL=1e-3, xTrain=X[train,], yTrain=Y[train],xTest=X[test,], yTest=Y[test])

foo <- stan(file="C:/Users/David Elkind/Desktop/RConvFn/Bayesian SVM/Bayes SVM - SMO.stan", data=iris_test, iter=200, chains=2, cores=7)

#############################################

C <- extract(foo, "C")[[1]]
gamma <- extract(foo, "gamma")[[1]]
beta1 <- extract(foo, "beta1")[[1]]
beta0 <- extract(foo, "beta0")[[1]]
z <- extract(foo, "z")[[1]]

bar <- cbind(ln_C=log(C), ln_gamma=log(gamma), beta0, beta1)
cor(bar)

print(foo, c("gamma", "C", "beta", "tau", "Omega", "mu", "lp__"))

pairs(foo, pars=c("beta", "beta", "C", "gamma", "lp__"))
pairs(foo, pars=c("gamma", "C"))
pairs(foo, pars=c("beta0", "rho"))
pairs(foo, pars=c("beta1", "C"))

pairs(foo, pars=c("alpha", "z"))
