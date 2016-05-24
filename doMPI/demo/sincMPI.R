library(doMPI)

# create and register a doMPI cluster if necessary
if (!identical(getDoParName(), 'doMPI')) {
  cl <- startMPIcluster(count=2)
  registerDoMPI(cl)
}

# define the grid over which to compute the sinc function
x <- seq(-20, 20, by=0.25)

# compute the sinc function in parallel
v <- foreach(y=x, .combine="cbind") %dopar% {
  r <- sqrt(x^2 + y^2) + .Machine$double.eps
  sin(r) / r
}

# display the results with a perspective plot
persp(x, x, v)
