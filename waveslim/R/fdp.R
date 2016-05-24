fdp.mle <- function(y, wf, J=log(length(y),2))
{
  fdpML <- function(d, y) {
    
    y.dwt <- y[[1]]
    n <- y[[2]]
    J <- y[[3]]

    ## Establish the limits of integration for the band-pass variances
    a <- c(1/2^c(1:J+1), 0)
    b <- 1/2^c(0:J+1)

    ## Define some useful parameters for computing the likelihood
    length.j <- n / c(2^(1:J), 2^J)
    scale.j <- c(2^(1:J+1), 2^(J+1))

    ## Initialize various parameters for computing the approximate ML
    bp.var <- numeric(J+1)

    ## Compute the band-pass variances according to d
    omega.diag <- NULL
    for(j in 1:(J+1)) {
      bp.var[j] <- integrate(fdp.sdf, a[j], b[j], d=d)$value
      omega.diag <- c(omega.diag, scale.j[j] * rep(bp.var[j], length.j[j]))
    }
    
    ## Compute approximate maximum likelihood 
    n * log(sum(y.dwt^2 / omega.diag) / n) +
      sum(length.j * log(scale.j * bp.var)) 
  }

  n <- length(y)
  y.dwt <- as.vector(unlist(dwt(y, wf, n.levels=J)))

  ## Compute MLE of d (limited to stationary region)
  result <- optimize(fdpML, interval=c(-0.5,0.5), maximum=FALSE,
                     y=list(y.dwt, n, J))

  ## Compute MLE of sigma_epsilon^2
  a <- c(1/2^c(1:J+1), 0)
  b <- 1/2^c(0:J+1)
  length.j <- n / c(2^(1:J), 2^J)
  scale.j <- c(2^(1:J+1), 2^(J+1))
  bp.var <- numeric(J+1)
  omega.diag <- NULL
  for(j in 1:(J+1)) {
    bp.var[j] <- integrate(fdp.sdf, a[j], b[j], d=result$minimum)$value
    omega.diag <- c(omega.diag, scale.j[j] * rep(bp.var[j], length.j[j]))
  }
  sigma2 <- sum(y.dwt^2 / omega.diag) / n

  list(parameters=c(result$minimum, sigma2), objective=result$objective)
}

