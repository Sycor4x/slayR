## This is a convenience function for computing the geodesics on a sphere with radius r by using the Vinchenty formula, which has the desirable property of being very precise in all use cases, including antipodal and collinear points, as well as having relatively low error for points without a high amount of decimal precision. (These are cases where other methods, like the law of haversines and spherical law of cosines, perform poorly.)
## The use case that I had in mind is to get a reasonable estimate of the distance between two points {x,y} on the Earth given by latitudes and longitudes in decimal degres. These will not be precise because (1) topography is ignored: hills, valleys, etc are all "smoothed over" and (2) this implementation assumes a spherical model: the Earth is not truely a sphere and (3) the default value of r is the mean radius of the Earth at the equator, so for points that are nearer to the poles, this will be slightly less accurate than for two points near the equator.
## Nonetheless, the total error in this approximation should be on the order of +/- 1%, which I think is resonable for the general use-cases in our practice, where the increased precision of more accurate methods will likely be swamped by statistical error elsewhere in the model.

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
