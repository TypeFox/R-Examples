dgammanc <- function(x, alpha, delta)
{
  n <- length(x) 
  if (n == length(alpha))   xx <- cbind(x, alpha) else
    if (length(alpha) == 1) xx <- cbind(x, rep(alpha, times = n))
  if (n == length(delta))   xx <- cbind(xx, delta) else
    if (length(delta) == 1) xx <- cbind(xx, rep(delta, times = n))
  dtched <- function(xx) return(dncgamma(xx[1], xx[2], xx[3]))
  d <- apply(xx, 1, dtched)
  return(d)
}

pgammanc <- function(x, alpha, delta)
{
  n <- length(x) 
  if (n == length(alpha))   xx <- cbind(x, alpha) else
    if (length(alpha) == 1) xx <- cbind(x, rep(alpha, times = n))
  if (n == length(delta))   xx <- cbind(xx, delta) else
    if (length(delta) == 1) xx <- cbind(xx, rep(delta, times = n))
  dtched <- function(xx) return(pncgamma(xx[1], xx[2], xx[3]))
  p <- apply(xx, 1, dtched)  
  return(p)
}

qgammanc <- function(p, alpha, delta)
{
  n <- length(p) 
  if (n == length(alpha))   xx <- cbind(p, alpha) else
    if (length(alpha) == 1) xx <- cbind(p, rep(alpha, times = n))
  if (n == length(delta))   xx <- cbind(xx, delta) else
    if (length(delta) == 1) xx <- cbind(xx, rep(delta, times = n))
  dtched <- function(xx) return(inxncgamma(xx[1], xx[2], xx[3]))
  q <- apply(xx, 1, dtched)  
  return(q)
}

deltagammanc <- function(x, alpha, p)
{
  n <- length(x) 
  if (n == length(alpha))   xx <- cbind(x, alpha) else
    if (length(alpha) == 1) xx <- cbind(x, rep(alpha, times = n))
  if (n == length(p))   xx <- cbind(xx, p) else
    if (length(p) == 1) xx <- cbind(xx, rep(p, times = n))
  dtched <- function(xx) return(indncgamma(xx[1], xx[2], xx[3]))
  delta  <- apply(xx, 1, dtched)  
  return(delta)
}