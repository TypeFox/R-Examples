`DDX` <-
function(x)
{
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  x <- rev(sort(x))
  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  #initial parameterization for nlm fx: estimates under yule model
  s1 <- yuleint2(x, x[1], 0)

  ddfunc <- function(r, v)
  {
    -(sum(log(2:(N-1))) + (N-2)*log(r) - sum(v*log(2:(N-1)))
    - sum(((2:N)^(1-v))*z *r) )
  }
  #res <- suppressWarnings(nlm(function(p) ddfunc(p[1], p[2]), c(s1$smax, 0), hessian = TRUE))
  res <- suppressWarnings(nlm(function(p) ddfunc(p[1], p[2]), c(s1[1], 0), hessian = TRUE))
  #may want to recode this to use 'optim' rather than 'nlm'
  aic <- 2*res$minimum + 4
  summ <- structure(list(LH = -res$minimum, aic = aic, r1 = res$estimate[1], xparam = res$estimate[2]));
  return(summ);

}

