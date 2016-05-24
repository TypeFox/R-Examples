`DDL` <-
function(x)
{
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  #calculates likelihoods under DD logistic model
  x <- rev(sort(x))
  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  ddfunc <- function(r, k)
  {
    -(sum(log(2:(N-1))) + (N-2)*log(r) + sum(log(1-((2:(N-1))/k)))
    - sum((2:N)*r*z) + sum(z*r*(2:N)^2)/k)
  }
   res <- suppressWarnings(nlm(function(p) ddfunc(p[1], p[2]), c(.5, N*2), hessian = TRUE))
   #may want to recode this to use 'optim' rather than 'nlm'
   aic <- 2*res$minimum + 4
   summ <- structure(list(LH = -res$minimum, aic = aic, r1 = res$estimate[1], kparam = res$estimate[2]))
   return(summ)
}

