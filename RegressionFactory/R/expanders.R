regfac.expand.1par <- function(beta, X, y, fbase1, fgh=2, ...) {
  # obtain base distribution derivatives
  ret <- fbase1(X%*%beta,y,fgh, ...)
  # expand base derivatives
  f <- sum(ret$f)
  if (fgh==0) return (f)
  g <- t(X)%*%ret$g
  if (fgh==1) return (list(f=f, g=g))
	xtw <- 0*X
	for (k in 1:ncol(X)) xtw[,k] <- X[,k]*ret$h # TODO: convert for loop to sapply
	h <- t(xtw)%*%X
  return (list(f=f, g=g, h=h))
}

regfac.expand.2par <- function(coeff, X, Z=matrix(1.0, nrow=nrow(X), ncol=1), y, fbase2, fgh=2, block.diag=FALSE, ...) {
  # extracting coefficients of X and Z
  K1 <- ncol(X); K2 <- ncol(Z)
  beta <- coeff[1:K1]
  gamma <- coeff[K1 + 1:K2]
  
  # obtain base distribution derivatives
  ret <- fbase2(X%*%beta, Z%*%gamma, y, fgh, ...)

  # expand base derivatives; TODO: vectorize for loops
  # function
  f <- sum(ret$f)
  if (fgh==0) return (f)
  # first derivative
  g <- c(t(X)%*%ret$g[,1],t(Z)%*%ret$g[,2])
  if (fgh==1) return (list(f=f, g=g))
  # second derivative
	h <- array(0, dim=c(K1+K2, K1+K2))
	# XX block
	xtw <- 0*X
	for (k in 1:K1) xtw[,k] <- X[,k]*ret$h[,1]
	h[1:K1, 1:K1] <- t(xtw)%*%X
	# ZZ block
	ztw <- 0*Z
	for (k in 1:K2) ztw[,k] <- Z[,k]*ret$h[,2]
	h[K1 + 1:K2, K1 + 1:K2] <- t(ztw)%*%Z
	# XZ and ZX blocks
	if (!block.diag) {
	  ztw2 <- 0*Z
	  for (k in 1:K2) ztw2[,k] <- Z[,k]*ret$h[,3]
	  h[K1 + 1:K2, 1:K1] <- t(ztw2)%*%X
	  h[1:K1, K1 + 1:K2] <- t(h[K1 + 1:K2, 1:K1])
	}
	
  return (list(f=f, g=g, h=h))
}




