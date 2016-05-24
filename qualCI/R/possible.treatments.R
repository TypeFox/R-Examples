possible.treatments <- function(d) {
  dat <- factor(d)
  N <- length(dat) # total number of units in set
  n <- tabulate(dat) # count of ctrl and treated units
  ng <- length(n) # number of values (should always be 2)
  if(ng!=2) stop("Need 2 treatment levels per set.")
  foo2 <- combn(N,n[2]) # these are essentially indices for control
  #out <- matrix(NA, nrow=N, ncol=ncol(foo2))
  out <- sapply(1:ncol(foo2),function(c){
  	vec <- rep(0,N)
  	vec[foo2[,c]] <- 1
  	vec
  	})
  out <- list(possibleTreat = t(out), prob = rep(1/nrow(t(out)),nrow(t(out))))
}