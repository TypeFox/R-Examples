####################################################################
#  Original R function written by Jaya Satagopan (2004)
#  Fortran backbone and modification by Venkatraman Seshan
####################################################################

twoStagePower <- function(n=NULL, Cost=NULL, m=5000, mu=0.045, mu.loc=0.5,
                          p=0.10, f=NULL, relcost=1, true.needed=1, rho=0,
                          rho0=0, nsim=2000) {
  if(is.null(n) & is.null(Cost)) stop("Specify exactly one of n or Cost")
  if(any(mu.loc <= 0 | mu.loc >= 1)) stop("locations of the true markers should be between 0 and 1")
# convert locations to indices
  mu.loc <- round(m*mu.loc)
  ntrue <- length(mu.loc)
  if(true.needed > ntrue) stop(paste("true.needed (= ", true.needed, ") exceeds number of true markers (= ", ntrue, ")", sep=""))
# calculate maximum and first stage sample size under Cost constraint
  if (is.null(n)) {
    if (is.null(f)) f <- 0.75
    n1 <- round(Cost * f)
    n <- round(Cost * (f + (1-f)/(relcost*p)))
    out <- c(Cost, p, f, 0)
    names(out) <- c("Cost", "p", "f", "power")
  } else {
# first stage sample size when maximum n is specified
    if (is.null(f)) f <- 0.5
    n1 <- round(n * f)
    out <- c(n, p, f, 0)
    names(out) <- c("n", "p", "f", "power")
  }
# other parameters
  n2 <- n - n1
  m1 <- round(m*p)
# calculate the mean vector
  mu.star <- rep(0,m)
  if (length(mu) == 1) {
    mu <- rep(mu, ntrue)
  } else if (length(mu) != ntrue) {
    stop("mu should be either a single number or of same length as mu.loc")
  }
  if(rho > 0) {
    tt <- matrix(0, nrow=ntrue, ncol=m)
    for(i in 1:ntrue)  tt[i,] <- mu[i] * rho^abs((1:m)-mu.loc[i])
    mu.star <- apply(tt,2,sum)
  } else {
    mu.star[mu.loc] <- mu
  }

# compute power
  zzz <- .Fortran("pwr2stg",
                  m=as.integer(m),
                  m1=as.integer(m1),
                  n1=as.integer(n1),
                  n2=as.integer(n2),
                  x1=double(m),
                  x2=double(m),
                  muvec=as.double(mu.star),
                  nmu=as.integer(ntrue),
                  muloc=as.integer(mu.loc),
                  rho0=as.double(rho0),
                  rho1=as.double(rho),
                  nsel=as.integer(true.needed),
                  nsim=as.integer(nsim),
                  pow=double(1),
                  tmp=double(m),
                  ord1=integer(m),
                  rnk1=integer(m))
  out[4] <- zzz$pow
  out
}
