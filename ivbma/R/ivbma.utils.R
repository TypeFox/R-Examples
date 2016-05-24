rinvwish.old <- function (v, S) 
{
	S <- solve (S)
	if (!is.matrix(S)) 
		S <- matrix(S)
	if (nrow(S) != ncol(S)) {
		stop(message = "S not square in rwish().\n")
	}
	if (v < nrow(S)) {
		stop(message = "v is less than the dimension of S in rwish().\n")
	}
	p <- nrow(S)
	CC <- chol(S)
	Z <- matrix(0, p, p)
	diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
	if (p > 1) {
		pseq <- 1:(p - 1)
		Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * (p - 1)/2)
	}
	ret <- crossprod(Z %*% CC)

	return(solve(ret))
}

rinvwish <- function (delta, D) 
{
  ##Bartlet Decomposition
  T <- chol(solve(D))
  p <- dim(D)[1]
  Psi <- matrix(0, p, p)
  for (i in 1:p) Psi[i, i] <- sqrt(rchisq(1, delta + p - i))
  if (p > 1)
    {
      for (i in 1:(p - 1))
        {
          for (j in (i + 1):p)
            {
              Psi[i, j] <- rnorm(1)
            }
        }
    }
  Psi <- Psi %*% T
  return(solve(t(Psi) %*% Psi))
}

rmvnorm <- function(mu, Sigma)
{
  p <- dim(Sigma)[2]
  z <- rnorm(p)
  x <- z %*% chol(Sigma)
  return(x + mu)
}

rmvnorm.precision <- function(mu,K)
  {
    p <- dim(K)[2]
    Q <- chol(K)
    z <- rnorm(p)
    x <- backsolve(Q,z) + mu
    return(x)
  }

mydet <- function(A)
  {
    return(2 * sum(log(diag(chol(A)))))
  }
