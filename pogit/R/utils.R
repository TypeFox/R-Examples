############################################################################## #
#### ------------- miscellaneous functions for pogit package -------------- ## #
############################################################################## #

# ---------------------------------------------------------------------------- #
#                         -- thins MCMC samples --
# ---------------------------------------------------------------------------- #
# used in: summary.pogit(), plot.pogit() 

thinMCMC <- function(x, start=NULL, mcmc=x$mcmc){
  if (is.null(x)) xth <- NULL
  if (!all(sapply(x, is.numeric))) stop("non-numeric values in 'x'")
  if (is.matrix(x)){
    if (is.null(start)) start <- 1
    idx <- seq(start, mcmc$nmc, mcmc$thin)
    xth <- x[idx, , drop = FALSE]
  }
  return(xth)
}

# ---------------------------------------------------------------------------- #
#   -- computes hightest posterior density (HPD) interval of MCMC samples --
# ---------------------------------------------------------------------------- #
# used in: summary.pogit()

hpdMCMC <- function(x, alpha = 0.05, ...){
  if (!is.vector(x)) x <- as.vector(x)
  nr <- NROW(x)
  l <- ceiling(alpha * nr)
  lb   <- sort(x)[1:l]
  ub   <- sort(x)[(nr - l + 1):nr]  
  minl <- min(ub - lb)
  lower <- lb[which((ub - lb)==minl)][1]
  upper <- ub[which((ub - lb)==minl)][1] 
  return(c(lower = lower, upper = upper))    
}


# ---------------------------------------------------------------------------- #
#    -- computes integrated autocorrelation times (IAT) of MCMC samples --
# ---------------------------------------------------------------------------- #
## used in: summary.pogit()

iatMCMC <- function(x){ 
  if (!is.vector(x)) x <- as.vector(x)
  nr   <- NROW(x)
  nlag <- max(3, floor(nr/2))  # maximum lag
  mu <- mean(x)
  vx <- var(x)
  c  <- matrix(0, nlag, 1) 
  
  # empirical autocovariances
  for (i in 1:nlag){
    x1 <- x[(i+1):nr] - mu
    x2 <- x[1:(nr-i)] - mu
    c[i] <- t(x1)%*%x2/(nr*vx)
  }
  
  index <- seq(2, nlag - 1, 2)
  # sums of adjacent pairs of autocovariances are positive
  Ga    <- c[index] + c[index + 1] 
  h1    <- which(Ga < 0)
  if (length(h1)==0){
    m1 <- (nlag - 1)/2    
  } else {
    m1 <- min(h1) - 1  
  }
  
  # first negative Gamma_m gives the initial positive sequence estimator
  # m is chosen to be the largest integer such that hat_Gamma>0
  if (m1 > 0){
    dGa <- diff(Ga)
    h2  <- which(dGa > 0)
    if (length(h2)==0){
      m <- m1
    } else {
      m <- min(m1, min(h2) - 1)
    }
  } else {
    m <- m1
  }

  iat <- 1 + 2*sum(c[1:(2*m + 1)])
  ess <- length(x)/iat
  return(c(IAT = iat, ESS = ess))
}



