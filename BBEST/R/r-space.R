##################################################################
#
#
# FUNCTIONS TO WORK WITH R-SPACE
#
# Authors: Charles R. Hogg III, Anton Gagin
#


#####################################################################################
# noise.cov.r(r1, r2, Q, sigma)
  # Computes the covariance between two points in G(r) due to noise in S(Q)
  # (the latter is assumed i.i.d. Gaussian).
  #
  # Args:
  #   r1:  (numeric) One r-value to consider
  #   r2:  (numeric) The other r-value to consider 
  #   Q:  (numeric vector) The Q-values where we have data.
  #   sigma:  (numeric) The standard deviation of the noise in S.
  #
  # Returns:
  #   (numeric) The covariance between G(r1) and G(r2) due to noise in S(Q).
noise.cov.r <- function(r1, r2, Q, sigma) {
  delta <- diff(Q)
  delta <- c(delta[1], delta)
  f.sum <- sum((2*sigma * Q * delta/pi ) ^ 2 * sin(Q * r1) * sin(Q * r2))  
  return (f.sum)
}

#####################################################################################
# noise.cov.vector.r(r1, r2, Q, sigma)
noise.cov.vector.r <- Vectorize(noise.cov.r, vectorize.args=c("r1", "r2"))

#####################################################################################
# noise.cov.matrix.r(r, Q, sigma)
#
  # Computes the covariance matrix in G(r) due to noise in S(Q) (the latter is
  # assumed i.i.d. Gaussian).
  #
  # Args:
  #   r:  (numeric vector) The r-values where we evaluate this covariance.
  #   Q:  (numeric vector) The Q-values where we have data.
  #   sigma:  (numeric) The standard deviation of the noise in S.
  #
  # Returns:
  #   A (N x N) matrix M, where N is the length of r, such that M[i, j] gives
  #   the covariance between r[i] and r[j].
noise.cov.matrix.r <- function(r, Q, sigma) {
  N <- length(r)
  M <- matrix(nrow=N, noise.cov.vector.r(r1=rep(r, each=N), r2=rep(r, N), Q=Q, sigma=sigma))
  return (M)
}

sineFT.matrix <- function(Q, r) {
  # Computes the matrix which converts a function at the given Q-points into
  # r-space using a Fourier sine transform (via Simpson equation).
  #
  # Args:
  #   Q:  (numeric vector) The Q-values where the function is evaluated.
  #   r:  (numeric vector) The r-values where the function is evaluated.
  #
  # Returns:
  #   A (N.r x N.Q) numeric matrix which effects the Q-to-r sine Fourier
  #   transform.
  dQ <- Dx(Q)
  N.Q <- length(Q)
  return(sapply(X=1:N.Q, Q=Q, dQ=dQ, r=r, 
           FUN=function(i, Q, dQ, r) {
             factor <- 2*(i%%2+1)/3
		     if(i==1 || i==length(Q))
		       factor <- 1/3	 
            2/pi * Q[i] * sin(Q[i]*r) * dQ[i]*factor
          }))
}

sineFT <- function(f.Q, Q, r) {
  # Computes the FT at the given Q-points into
  # r-space using a Fourier sine transform.
  N.Q <- length(Q)
  N.r <- length(r)
  dQ <- Dx(Q)
  f.r <- 0 * r
  
  for (i in seq(2, N.Q-1, 2)) {
    f.r <- f.r + 2/pi *Q[i]*sin(Q[i] * r) * f.Q[i] * dQ[i] * 4/3
  }
  for (i in seq(3, N.Q-1, 2)) {
    f.r <- f.r + 2/pi *Q[i]*sin(Q[i] * r) * f.Q[i] * dQ[i] * 2/3
  }
  
  f.r <- f.r + 2/pi*(Q[1]*sin(Q[1]*r)*f.Q[1]*dQ[1] + Q[N.Q]*sin(Q[N.Q]*r)*f.Q[N.Q]*dQ[N.Q] )*1/3
	
  return (f.r)
}


invert.order <- function(i) {
  # Inverts the result of the 'order' function.
  #
  # Args:
  #   i:  Numeric vector with a permutation of the integers from 1:length(i).
  #
  # Returns:
  #   Numeric vector r such that r[i] == 1:length(i).
  r <- i + NA
  for (j in 1:length(i)) {
    r[i[j]] <- j
  }
  return (r)
}

Dx <- function(x) {
  # Compute the width for each x (specifically, its Voronoi cell size) to aid
  # in numerical integration.
  #
  # Args:
  #   x:  A sorted numeric vector of x-values
  n <- length(x)
  i <- order(x)
  x.sort <- x[i]
  dx <- diff(c(x.sort[1], 0.5 * (x.sort[-1] + x.sort[-n]), x.sort[n]))
  return (dx[invert.order(i)])
}
