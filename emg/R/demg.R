demg <- function(x, mu=0, sigma=1, lambda=1, log=FALSE)
{
  l <- max(length(x), length(mu), length(sigma), length(lambda))
  
  x     <- rep(x,     times=ceiling(l/length(x)),     length.out=l)
  mu    <- rep(mu,    times=ceiling(l/length(mu)),    length.out=l)
  sigma <- rep(sigma, times=ceiling(l/length(sigma)), length.out=l)
  lambda<- rep(lambda,times=ceiling(l/length(lambda)),length.out=l)
  
  if(min(sigma) <= 0.0)       {stop("Sigma must be greater than zero") }
  if(min(lambda) <= 0.0)      {stop("Lambda must be greater than zero")}

  erfc <- pnorm((mu+lambda*sigma*sigma-x)/ sigma , lower.tail = FALSE, log.p=log)
  
  if(log)
  {
    result <- lambda/2 * (2*mu+lambda*sigma*sigma-2*x) + Re(erfc) + log(lambda)
  } else {
    result <- exp(lambda/2 * (2*mu+lambda*sigma*sigma-2*x)) * Re(erfc) * lambda
  }
  
  # Too small to compute, zero it
  result[is.nan(result)] <- 0

  result
}

