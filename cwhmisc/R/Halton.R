HS247 <- function(K,N,R,P=rep(0,K)) {
  # Halton, J.H., Smith, G.G., 1961. Algorithm 247
  # Radical-inverse quasi-random point sequence
  #  Computes a sequence of N quasi-random points 
  #  lying in the K-dimensional unit cube given by
  #  0 < x_i < 1, i = 1,2,..,K. The i-th component
  #  of the m-th point is stored in Q[m,i].
  #  The sequence is initiated by a zero-th point 
  #  stored in P, and each component sequence is
  #  iteratively generated with parameter R[i].
  #  E is a positive error-parameter.
  #  K, N, E, P[i], R[i], i=1..K, are to be given.
  Q <- matrix(nrow=N,ncol=K)
  for (ii in 1:K) {
    r <- 1.0/R[ii]
    for (nn in 1:N) {
      f <- if (nn > 1) 1.0 - Q[nn-1,ii] else 1.0 - P[ii]
      g <- 1.0
      h <- r
      while (f - h < 0.1) { # E = 0.1, nuisance parameter
        g <- h
        h <- h*r
      }
      Q[nn,ii] <- g + h - f
    }
  }
  return(Q)
}
