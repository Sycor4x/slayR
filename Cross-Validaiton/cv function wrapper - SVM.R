##	This function is designed for use inside of a CV routine in which alternative hyperparameter tuples will be compared.
##	Step 1: linear SVM, using e1071
library(cvTools)
inv_logit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x)-log(1-x)
inv_xform <- function(x) 1-inv_logit(x)

trapArea <- function(x,y){
	if(length(y)!=length(x)) stop("User error: length(y)!=length(x).")
	n <- length(y)
	width <- diff(x)
	height <- (y[1:(n-1)]+y[2:n])/2
	area <- sum(width*height)
	return(area)
}

AUC <- function(preds, label){
	if(length(preds)!=length(label)) stop("Error: length(preds)!=length(label).")
	if( !all(label==1|label==-1) ) stop("Error: label may only contain 0 or 1.")
	label <- ifelse(label==-1,0,1)
	o <- order(-preds)
	N <- length(preds)
	preds <- preds[o]
	label <- label[o]
	TP <- integer(N)
	FP <- integer(N)
	TP[1] <- xor(FALSE,	label[1])
	FP[1] <- xor(TRUE,	label[1])
	for(n in 2:N){
		TP[n] <- TP[n-1]+xor(FALSE, label[n])
		FP[n] <- FP[n-1]+xor(TRUE,	label[n])
	}
	FN <- sum(label)-TP
	TN <- N-sum(label)-FP
	TPR <- c(0,TP/(TP+FN))
	FPR <- c(0,FP/(FP+TN))
	AUC <- trapArea(x=FPR, y=TPR)
	return(AUC)
}

svm2 <- function(pars, X_train, Y_train){
	if( !all(c("ln_gamma", "ln_C")%in%names(pars) ) ) stop("User error: must supply named element 'ln_C' to function svm2 in argument pars. If this is called inside of slayR, check that you've declared rownames for argument 'box'.")
	library(e1071)
	C <- exp(pars[names(pars)=="ln_C"])
	gamma <- exp(pars[names(pars)=="ln_gamma"])
	svmOut <- svm(x=X_train, y=as.factor(Y_train), scale=FALSE, type="C-classification", kernel="radial", gamma=gamma, cost=C, cache=1024)
	return(svmOut)
}

applyKernel <- function(newX, FUN, ...){
  foreach(i=1:nrow(newX))%dopar%FUN(newX[i,], ...)
}

predict2 <- function(svmObj, X_test, Y){
	predOut <- predict(object=svmObj, newdata=X_test, decision.values=TRUE)
	decision <- attr(predOut, "decision.values")
	preds <- decision*Y[1]	## This vexing e1071 bug flips the sign of the decision function on the basis of whatever observation comes first
	return(preds)
}

cvSVM <- function(whatFold, X_in, Y, folds, pars){	
	r <- whatFold[1]
	k <- whatFold[2]
	# train <- folds$subsets[folds$which!=k,r] 
	test <- folds$subsets[folds$which==k,r]
	
	X_train <- X_in[-test,];			Y_train <- Y[-test];
	X_test <- X_in[test,];				Y_test <- Y[test]

	svmOut <- svm2(pars=pars, X_train=X_train, Y_train=Y_train)
	preds <- predict2(svmObj=svmOut, X=X_test, Y=Y_train)
	perf <- AUC(preds=preds, label=Y_test)
	return(perf)
}

fitSVM <- function(pars, Xdata, Y, folds){
	r <- folds$R
	k <- folds$K
	cvGrid <- as.matrix(expand.grid(r=1:r, k=1:k))
	cvAUC <- apply(cvGrid, MARGIN=1, FUN=cvSVM, X_in=Xdata, folds=folds, pars=pars, Y=Y)
	cvAUC <- sapply(cvAUC, function(x) median(c(x,1-1e-6,1e-6)) )
	logit_cvAUC <- logit(1-cvAUC) ## take negative because we want to maximize AUC, but we're using a minimization algorithm
	print(mean(logit_cvAUC))
	return(data.frame(y_mean=mean(logit_cvAUC), y_var=var(logit_cvAUC)))
}

library(data.table)
X <- data.table(iris)
X[,Species:=ifelse(X[,Species]=="virginica", 1, -1)]
X <- rbind(X, X[Species==1]) ## double the size of the minority class so that we have a 50/50 class balance. This way CV will probably have all members of the same class
Y <- X[,Species]
X <- as.matrix(X[,list(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)])
X <- scale(X)

foo <- fitSVM(pars=c(ln_C=5, ln_gamma=-2), X=X, Y=Y, folds=cvFolds(nrow(X), R=10, K=5), svmKernel=GaussianSVMKernel)

svmBox <- data.frame(min=c(log(2^-5),log(2^-15)),max=c(log(2^15),log(2^3)))
rownames(svmBox) <- c("ln_C", "ln_gamma")
bar <- slayR(f=fitSVM, box=svmBox, Xdata=X, Y=Y, folds=cvFolds(nrow(X), R=2, K=5), svmKernel=GaussianSVMKernel, kernel="Matern", iter=20)

plot(bar$xfx_init)


