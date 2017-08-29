library(rstan)
library(MASS)

library(sp)
library(rworldmap)
library(dplyr)
library(rgeos)

# Outlines for the US and Mexico
`%\\%` <- setdiff

US <- countriesLow %>% subset(NAME == "United States")
Mexico <- countriesLow %>% subset(NAME == "Mexico")

# Puff up Mexico by 25 miles, and then take only the overlap with the US
mexicoBuffer <- gBuffer(Mexico, width = .36)
borderBuffer <- gIntersection(mexicoBuffer, US)

# Make 1000+/- hexagonal grid cells
hexgrid <- spsample(borderBuffer, 1000, "hexagonal")

# Take a random sample of them
X <- as.matrix(sample_n(as.data.frame(hexgrid@coords), size = 100))
plot(X, xlab="Latitude", ylab="Longitude")

o <- order(X[,1])
X <- X[o,]

# This program is designed for inference on this model:
# Observations:
# X, the n lat/long locations of the observations
# Y, the counts of events at the lat/long locations
# Model:
# Y[i] ~ Poisson(z[i])
# ln z ~ Normal(mu*1[n], Sigma) This is a multivariate normal distribution with a mean vector constant at mu for all n entries, and Sigma a covariance matrix defined below
# Sigma[i,j] = Matern_5/2(geodesic(X[i,], X[j,]))
# geodesic is the (shortest) distance between two points on the surface of a sphere

inv_logit <- function(x) 1/(1+exp(-x))

geodesic <- function(x, y, r=6378.1370){
  phi <- c(x[1], y[1])*pi/180
  lambda <- c(x[2], y[2])*pi/180
  delta_lambda <- abs(diff(lambda))

  term1 <- (cos(phi[2])*sin(delta_lambda))^2
  term2 <- (cos(phi[1])*sin(phi[2])-sin(phi[1])*cos(phi[2])*cos(delta_lambda))^2
  denominator <- sin(phi[1])*sin(phi[2])+cos(phi[1])*cos(phi[2])*cos(delta_lambda)

  delta_sigma <- atan2(sqrt(term1+term2), denominator)
  distance_out <- r*delta_sigma;
  return(distance_out)  
}

Matern_geo <- function(X, inv_rho, sigma_sq=0){
  K <- matrix(0, nrow(X), nrow(X))
  coef <- c(1, sqrt(5), 5/3)
  fact2 <- numeric(3)
  for(i in 1:(nrow(X)-1)){
    for(j in (i+1):nrow(X)){
      dist_rho_quotient <- inv_rho*geodesic(X[i,], X[j,])
      fact2 <- dist_rho_quotient^(0:2)
      K[i,j] <- exp(-sqrt(5)*dist_rho_quotient)*(coef%*%fact2)
      K[j,i] <- K[i,j]
    }
  }
  diag(K) <- 1+sigma_sq
  return(K)
}

K <- Matern_geo(X, inv_rho=0.01, sigma_sq=1e-6)  ## K+sigma_sq*I
all(eigen(K)$values>0)  ## should be true -- if not, the culprit is plausibly numerical singularity. Increase sigma.

z <- mvrnorm(n=1, mu=rep(1.5, nrow(K)), Sigma=K)
y <- rpois(length(z), lambda=exp(z))

plot(X, xlab="Latitude", ylab="Longitude", type="n")
text(X[,1], X[,2], y)

N <- 80
M <- nrow(X)-N
train <- sample(1:(M+N), size=N, replace=FALSE)
test <- (1:(M+N))%\%train
x_train <- X[train,]
y_train <- y[train]
x_test <- X[test,]
y_test <- y[test]
stanTest <- list(M=M, N=N, x_train=x_train, x_test=x_test, y_train=y_train, y_test=y_test)

all(train%in%test) ## should be disjoint

bar <- stan(file="C:/Users/David Elkind/Desktop/RConvFn/Geospatial Modeling/Spatial Matern model.stan", data=stanTest, cores=4, iter=2000, chains=4)

# foo <- stan_model(file="C:/Users/David Elkind/Desktop/RConvFn/Geospatial Predictive Modeling with Kernels/Spatial Matern model.stan")
# bar <- optimizing(foo, as_vector=FALSE)

print(bar, c("inv_rho", "mu", "sigma_sq", "rmse"))

z_train <- apply(extract(bar, "z_train")[[1]], 2, mean)
plot(z[train], z_train, xlim=c(-5,5), ylim=c(-5,5))
z_diagnostic <- lm(z_train~z[train])
abline(a=0,b=1)
abline(z_diagnostic)

z_test <- colMeans(extract(bar, "z_test")[[1]])
plot(z[test], z_test, xlim=c(-5,5), ylim=c(-5,5))
z_diagnostic <- lm(z_test~z[test])
abline(a=0,b=1)
abline(z_diagnostic)

plot(density(extract(bar, "rmse")[[1]], from=0))
abline(h=mean(y_test))