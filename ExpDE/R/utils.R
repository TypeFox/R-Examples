# Generate an matrix of uniformly distributed values in the interval (0,1) with
# the same dimension as M
randM <- function(M) {
  matrix(stats::runif(prod(dim(M))), 
         ncol = ncol(M))
}

# Denormalize population
denormalize_population <- function(probpars, Pop){
  # Denormalize population
  LL <- matrix(rep(probpars$xmin, nrow(Pop)),
               ncol = ncol(Pop),
               byrow = TRUE)
  UL <- matrix(rep(probpars$xmax, nrow(Pop)),
               ncol = ncol(Pop),
               byrow = TRUE)
  return(LL + Pop*(UL-LL))
}