# stationary mean function
stationary <- function(t) { cbind( array(1,length(t)) ) }

# Fourier series of period T and degree k
# Can we do this with an FFT? For irregular data?
periodic <- function(t,k=1,P=24*60^2)
{
  u <- (2*pi/P)*t
  u <- u %o% (1:k)
  
  # rifle cosine and sine matrices
  u <- cbind( cos(u) , sin(u) )
  u <- u[,order(sequence(c(k,k)))]

  # add constant term
  u <- cbind( stationary(t) , u )
  
  return(u)
}

# polynomial
polynomial <- function(t,k=1)
{
  u <- sapply(t, function(s) { s^(1:k) } )
  u <- t(u)
  u <- cbind( stationary(t) , u )
  
  return(u)
}

# evaluate a mean function over a range of times
call.mean <- function(CTMM,t)
{
  par <- c(list(t=t),attr(CTMM$mean,"par"))
  do.call(CTMM$mean,par)
}

# calculate deterministic ms speed and variance
mspeed <- function(CTMM)
{
  if(CTMM$mean=="stationary")
  {
    MS <- 0
    VAR <- 0
  }
  else if(CTMM$mean=="periodic")
  {
    STUFF <- prepare.periodic(CTMM)
    omega <- STUFF$omega
    A <- STUFF$A
    COV <- STUFF$COV
    
    MS <- sum(omega^2*A^2)/2
    grad <- omega^2*A
    VAR <- grad %*% COV %*% grad
  }
  
  return(list(MS=MS,VAR=VAR))
}

# extract some useful info from periodic mean function parameters
prepare.periodic <- function(CTMM)
{
  P <- attr(CTMM$mean,"par")$P
  omega <- 2*pi/P
  
  # amplitudes and covariance
  A <- CTMM$mu
  COV <- CTMM$COV.mu
  k <- (nrow(A)-1)/2
  
  # harmonic numbers
  K <- c( 1:k , 1:k )
  K <- K[order(sequence(c(k,k)))]
  K <- c(0,K)
  
  # flatten block-vectors and block-matrices
  A <- array(A,2*(2*k+1))
  COV <- array(COV,c(2,2*k+1,2*k+1,2))
  COV <- aperm(COV,c(2,1,3,4))
  COV <- array(COV,c(2*(2*k+1),2*(2*k+1)))
  
  omega <- K*omega
  
  return(list(A=A,COV=COV,omega=omega))
}

# SVF of mean function
svf.mean <- function(CTMM)
{
  if(CTMM$mean=="stationary")
  {
    svf <- function(t) { 0 }
    VAR <- function(t) { 0 }
  }
  else if(CTMM$mean=="periodic")
  {
    STUFF <- prepare.periodic(CTMM)
    omega <- STUFF$omega
    A <- STUFF$A
    COV <- STUFF$COV
    
    svf <- function(t) { sum( 1/4 * A^2 * (1-cos(omega*t) )) }
    VAR <- function(t)
    {
      grad <- 1/2 * A * (1-cos(omega*t))
      return(grad %*% COV %*% grad)
    }
  }
  
  return(list(svf=svf,VAR=VAR))
}