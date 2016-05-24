rtruncdist <- function(n, ..., dist='norm', truncmin=-Inf, 
                       truncmax=Inf){  
  qtruncdist(runif(n), ..., dist=dist, truncmin=truncmin, 
             truncmax=truncmax)
} 