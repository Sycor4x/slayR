######
library(cvTools)
##	This function is designed for use inside of a CV routine in which alternative hyperparameter tuples will be compared.
##	Step 1: linear SVM, using e1071
library(cvTools)
inv_logit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x)-log(1-x)
inv_xform <- function(x) 1-inv_logit(x)

applyKernel <- function(newX, FUN, ...){
  foreach(i=1:nrow(newX))%dopar%FUN(newX[i,], ...)
}

fitSVM <- function(pars, Xdata, Y, folds){
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
		label <- ifelse(label==-1,0,1)
		if( !all(label==1|label==0) ) stop("Error: label may only contain 0 or 1.")
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

# X <- data.table(read.csv(file="C:/Users/David Elkind/Desktop/RConvFn/expanded iris data/newIris.csv", stringsAsFactors=FALSE))
# X <- X[sample(nrow(X), size=500)]
X <- data.table(iris)
X[,Species:=ifelse(X[,Species]=="versicolor", 1, -1)]
#X <- rbind(X, X[Species==1,])
Y <- X[,Species]
X <- as.matrix(X[,Species:=NULL])
X <- as.matrix(scale(X))

plot(data.frame(X), col=ifelse(Y<0, "blue", "red"))

windows()
plot(iris)

ln_C <- log(2^(seq(-5, 15, len=40+1)))
ln_gamma <- log(2^seq(-15, 3, len=40))
grid <- as.matrix(expand.grid(ln_C=ln_C, ln_gamma=ln_gamma))
folds <- cvFolds(nrow(X), R=25, K=5)

library(doParallel)
library(foreach)
cl <- makeCluster(7)
registerDoParallel(cl)
foo <- applyKernel(grid, fitSVM, folds=folds, Xdata=X, Y=Y)
stopCluster(cl)

foo <- do.call(rbind, foo)
bar <- cbind(grid, foo)
asdf <- unstack(bar, y_mean~ln_gamma)
library(fields)
image.plot(z=t(as.matrix(inv_xform(asdf))), x=ln_gamma, y=ln_C, ylab=expression(log(C)), xlab=expression(log(gamma)), main="AUC")

library(e1071)
svmTest1 <- svm(x=X, y=as.factor(Y), kernel="radial", type="C-classification", scale=FALSE)
pred1 <- attr(predict(svmTest1, newdata=X, decision.values=TRUE), "decision.values")
svmTest2 <- svm(x=X, y=as.factor(-Y), kernel="radial", type="C-classification", scale=FALSE)
pred2 <- attr(predict(svmTest2, newdata=X, decision.values=TRUE), "decision.values")
plot(pred1,pred2, main("The sign of Y[1] doesn't matter."))
abline(a=0,b=1, col='red')

##############################################

svmBox <- data.frame(min=c(log(2^-5),log(2^-15)),max=c(log(2^15),log(2^3)))
rownames(svmBox) <- c("ln_C", "ln_gamma")
bar <- slayR(f=fitSVM, box=svmBox, Xdata=X, Y=Y, folds=cvFolds(nrow(X), R=2, K=5), kernel="Matern", iter=10, init=10)

image.plot(z=t(as.matrix(inv_xform(asdf))), x=ln_gamma, y=ln_C, ylab=expression(log(C)), xlab=expression(log(gamma)), main="AUC", zlim=c(0.5,1))
points(bar$opt_path, col="black")
lines(bar$opt_path, col="black")
points(bar$xfx_init, col="grey")
points(bar$best, col="white", lwd=2)
