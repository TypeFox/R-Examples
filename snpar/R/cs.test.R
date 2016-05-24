cs.test <-
function(x, alternative = c("two.sided","increasing","decreasing"),
                    exact = TRUE, correct = TRUE)
{
  DNAME <- deparse(substitute(x))
  METHOD <- ifelse(exact,"Exact Cox-Stuart trend test",
                   "Approximate Cox-Stuart trend test")
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  x <- x[is.finite(x)]  
  n <- as.integer(length(x))
  if (n < 4L) 
    stop("not enough 'x' data") 
  c <- as.integer(ifelse(n%%2 == 0, n/2, (n - 1)/2))
  DIFF <- if (n%%2 == 0) x[1:c] - x[(c+1):n] else x[1:c] - x[(c+2):n]
  DIFF <- DIFF[DIFF != 0]
  c <- length(DIFF)
  s.p <- sum(DIFF > 0); s.m <- sum(DIFF < 0)
  alternative <- match.arg(alternative)
  if (correct)
  {
    if (s.m < length(DIFF)/2) 
      z = (s.m + 0.5 - length(DIFF)/2)/sqrt(length(DIFF)/4)
    else 
      z = (s.m - 0.5 - length(DIFF)/2)/sqrt(length(DIFF)/4)
  }
  else
    z = (s.m - length(DIFF)/2)/sqrt(length(DIFF)/4)
  p.val.mono <- ifelse(exact, 2*min(1 - pbinom(min(s.m,s.p),c,0.5),
                       pbinom(min(s.m,s.p),c,0.5)), 2*min(pnorm(z), 1 - pnorm(z)))
  p.val.incr <- ifelse(exact, pbinom(s.p,c,0.5), 1 - pnorm(z))
  p.val.decr <- ifelse(exact, pbinom(s.m,c,0.5), pnorm(z))
  p.val <- switch(alternative, two.sided = p.val.mono, increasing = p.val.incr, 
                  decreasing = p.val.decr)
  STAT <- switch(alternative, two.sided = min(s.m,s.p), increasing = s.p,
                 decreasing = s.m)
  STATISTIC <- setNames(STAT, switch(alternative, two.sided = "S", increasing = "S+",
                                     decreasing = "S-"))
  ALTERNATIVE <- switch(alternative,two.sided = "data have a monotonic trend", 
                        increasing = "data have an increasing trend",
                        decreasing = "data have a decreasing trend")
  RVAL <- list(data.name = DNAME, method = METHOD, alternative = ALTERNATIVE,
               p.value = p.val, statistic = STATISTIC)
  class(RVAL) <- "htest"
  return(RVAL)
}
