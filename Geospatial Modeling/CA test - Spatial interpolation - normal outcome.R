library(rstan)
library(MASS)
# Outlines for the US and Mexico
`%\\%` <- setdiff

library(sp)
library(dplyr)
library(rgeos)
library(rworldmap)

norm01 <- function(x, x_min, x_max){
  if(missing(x_min)&missing(x_max)) out <- (x-min(x))/(max(x)-min(x))
  if(!missing(x_min)&!missing(x_max)) out <- (x-x_min)/(x_max-x_min)
  if((missing(x_min)|missing(x_max))&!(missing(x_min)&missing(x_max))) stop("User error: Must supply both or neither of x_min, x_max.")
  return(out)
}

# Download this file: http://biogeo.ucdavis.edu/data/gadm2.8/rds/USA_adm1.rds
USStates <- readRDS("C:/Users/David Elkind/Desktop/RConvFn/Geospatial Modeling/Mexico Border/USA_adm1.rds")
CA <- USStates %>% subset(NAME_1 == "California")

# Drawn free-hand using http://dev.openlayers.org/examples/vector-formats.html
# http://code.boyles.cc/DrawPasteAndShow/
boundary <- readWKT("POLYGON((-115.54149627686 31.727828979492, -122.52880096436 37.748336791992, -124.28661346436 41.923141479492, -128.68114471436 41.945114135742, -125.78075408936 33.749313354492, -122.19921112061 30.783004760742, -117.98046112061 30.299606323242, -116.53026580811 30.761032104492, -115.54149627686 31.727828979492))")

# Puff up Mexico by 25 miles, and then take only the overlap with the US
borderBuffer <- CA %>%
  gBuffer(width = 1) %>%
  gDifference(countriesLow) %>%
  gIntersection(boundary)

# Make 1000 hexagonal grid cells
hexgrid <- spsample(borderBuffer, 1000, "hexagonal") %>%
  HexPoints2SpatialPolygons()

plot(hexgrid, bg="deepskyblue")
plot(countriesLow, col="gray", add=TRUE)
plot(USStates, add=TRUE)

# Write it to CSV
write.csv(coordinates(hexgrid), "C:/Users/David Elkind/Desktop/RConvFn/Geospatial Modeling/California Coastline/MPAT grid.csv", row.names = FALSE)

X <- coordinates(hexgrid)

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

Matern_geo <- function(X, inv_rho, eta_sq=1, sigma_sq=0){
  K <- matrix(0, nrow(X), nrow(X))
  coef <- c(1, sqrt(5), 5/3)
  fact2 <- numeric(3)
  for(i in 1:(nrow(X)-1)){
    for(j in (i+1):nrow(X)){
      dist_rho_quotient <- inv_rho*geodesic(X[i,], X[j,])
      fact2 <- dist_rho_quotient^(0:2)
      K[i,j] <- eta_sq*exp(-sqrt(5)*dist_rho_quotient)*(coef%*%fact2)
      K[j,i] <- K[i,j]
    }
  }
  diag(K) <- eta_sq+sigma_sq
  return(K)
}

K <- Matern_geo(X, inv_rho=0.001, eta_sq=2, sigma_sq=1e-6)  ## K+sigma_sq*I
all(eigen(K)$values>0)  ## should be true -- if not, the culprit is plausibly numerical singularity. Increase sigma.

y <- mvrnorm(n=1, mu=rep(1.5, nrow(K)), Sigma=K)

plot(X, xlab="Latitude", ylab="Longitude", type="n")
text(X[,1], X[,2], round(y))

N <- 9
M <- nrow(X)-N

A <- cbind(X, 1:nrow(X))
ndx <- A[X[,2]>32.5 & X[,2]<34.5,][,3]

train <- sort(sample(ndx, size=N, replace=FALSE))
test <- (1:(M+N))%\%train
x_train <- X[train,]
y_train <- y[train]
x_test <- X[test,]
y_test <- y[test]
stanTest <- list(M=M, N=N, x_train=x_train, x_test=x_test, y_train=y_train, y_test=y_test)

plot(hexgrid, col=grey(norm01(y, x_min=min(y), x_max=max(y))), border="transparent", main="True", bg="deepskyblue")
plot(hexgrid[train], col="red", add=TRUE, border="transparent")
plot(countriesLow, col="gray", add=TRUE)
plot(USStates, add=TRUE)

all(train%in%test) ## should be disjoint, i.e. FALSE

bar <- stan(file="C:/Users/David Elkind/Desktop/RConvFn/Geospatial Modeling/Spatial interpolation - normal outcome.stan", data=stanTest, cores=4, iter=2000, chains=4)

# foo <- stan_model(file="C:/Users/David Elkind/Desktop/RConvFn/Geospatial Predictive Modeling with Kernels/Spatial Matern model.stan")
# bar <- optimizing(foo, as_vector=FALSE)

print(bar, c("inv_rho", "mu", "sigma_sq", "rmse"))

y_hat <- colMeans(extract(bar, "y_hat")[[1]])
plot(y_test, y_hat)
y_diagnostic <- lm(y_test~y_hat)
abline(a=0,b=1)
abline(y_diagnostic)

centroid <- colMeans(x_train)
far <- apply(x_test, 1, geodesic, y=centroid)
plot(far, (y_test-y_hat)^2, xlab="Distance from Centroid of Sensors", ylab="Squared error")

resid <- -(y_test-y_hat)^2
resid_norm <- norm01(resid)
resid_proj <- rep(median(resid_norm), N+M)
resid_proj[test] <- resid_norm
plot(hexgrid, col=grey(resid_proj), border="transparent", main="Residuals", sub="white = best, black=worst", bg="skyblue")
plot(hexgrid[train], col="red", add=TRUE, border="transparent")
plot(countriesLow, col="gray", add=TRUE)
plot(CA, col="wheat1", add=TRUE)

y_hat_proj <- rep(median(y_hat), N+M)
y_hat_proj[test] <- y_hat
par(mfrow=c(1,2))
plot(hexgrid, col=grey(norm01(y, x_min=min(c(y, y_hat)), x_max=max(c(y,y_hat)))), border="transparent", main="True", bg="skyblue")
plot(hexgrid[train], col="red", add=TRUE, border="transparent")
plot(countriesLow, col="gray", add=TRUE)
plot(CA, col="wheat1", add=TRUE)
plot(hexgrid, col=grey(norm01(y_hat_proj, x_min=min(c(y, y_hat)), x_max=max(c(y,y_hat)))), border="transparent", main="Predicted", bg="skyblue")
plot(hexgrid[train], col="red", add=TRUE, border="transparent")
plot(countriesLow, col="gray", add=TRUE)
plot(CA, col="wheat1", add=TRUE)

plot(density(extract(bar, "rmse")[[1]], from=0))
abline(v=mean(y_test))