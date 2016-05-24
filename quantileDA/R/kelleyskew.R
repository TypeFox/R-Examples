kelleyskew<-function(x){
  #
  #  Compute Kelley's skewness measure for x
  #  NOTE: this procedure assumes no x values are missing
  #
  quarts <- as.numeric(quantile(x, probs = c(0.1, 0.5, 0.9)))
  num <- quarts[1] + quarts[3] - 2*quarts[2]
  denom <- quarts[3] - quarts[1]
  gskew <- num/denom
  gskew
}