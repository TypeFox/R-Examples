

rvdens <- function(n=1, FUN, range, unitprecision=10, ...) {
  # NAME
  #   rvdensity - Sample from a given univariate density using a grid approximation
  # ARGUMENTS
  #   n : number of independent random vector components to draw
  #   FUN : density function, must be vectorized
  #   range : range for the grid
  #   unitprecision : number of points per unit
  #   ... : other arguments passed to [FUN].
  #   
  grid <- seq(from=range[1], to=range[2], by=1/unitprecision)
  prob <- FUN(grid, ...)
  n.sims <- getnsims()
  s <- sample(grid, size=n*n.sims, prob=prob, replace=TRUE)
  noise <- runif(n*n.sims, -0.5/unitprecision, 0.5/unitprecision)
  rvsims(matrix(s+noise, nrow=n.sims))
}


