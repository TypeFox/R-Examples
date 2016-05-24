BBSim <- function(start, end, n, Sigma=diag(2), T, t0=0,seed=1){
  p <-   dim(Sigma)[1]
  start <- matrix(start,ncol=p)
  end   <- matrix(end,ncol=p)
#  if(seed>0){armadillo_set_seed(seed)}
#  if(seed<=0){armadillo_set_seed_random()}


  A <- .Call("BBSim",
             start,
             end,
             n,
             Sigma,
             T,
             t0,
             seed,
             p,
#             DUP=FALSE,
             PACKAGE = "BIPOD" )
  return(A)
}

