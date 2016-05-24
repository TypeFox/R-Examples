ddpareto <- function(x, theta, log = FALSE)
{
  if (theta <= 0 | theta>=1) {
    stop(paste("theta must be between 0 and 1!", "\n"))
  }
  n <- which(x < 0)
  x <- pmax(x, 0)
  p <- theta^(log(1+x))-theta^(log(2+x))
  if(length(n)>0) p[n] <- 0
  if(log) p <- log(p)
  p[is.nan(p)] <- 0
  if(!log) p <- pmin(pmax(p, 0), 1)
  p
}

pdpareto <- function(q, theta, lower.tail = TRUE, log.p = FALSE)
{
  if (theta <= 0 | theta>=1) {
    stop(paste("theta must be between 0 and 1!", "\n"))
  }
  ind <- (q < 0)
  q <- floor(q)
  temp <- 1 - theta^(log(2 + q))
  if(lower.tail==FALSE) temp <- 1 - temp
  if(any(ind)) temp[ind] <- 0 + 1*!lower.tail
  if(log.p) temp <- log(temp)
  if(!log.p) temp <- pmin(pmax(temp, 0), 1)
  temp 
}

qdpareto <- function(p, theta, lower.tail = TRUE, log.p = FALSE)
{
  if (theta <= 0 | theta>=1) {
    stop(paste("theta must be between 0 and 1!", "\n"))
  }
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  all.p <- pmax(floor(exp(log(1-p)/log(theta))-2),0)
  all.p[p==1] <- Inf
  all.p[p==0] <- 0
  all.p[(p>1)|(p<0)] <- NaN
  if(any(is.nan(all.p))) warning("NaNs produced")
  all.p
}

rdpareto <- function(n, theta)
{
  if (theta <= 0 | theta>=1) {
    stop(paste("theta must be between 0 and 1!", "\n"))
  }
  out <- floor(exp(rexp(n,-log(theta)))-1)
  out
}
