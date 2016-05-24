qn <- function(x)  
{ 
  if (NCOL(x) > 1)
    stop("'x' must be a numeric or univariate time series")
  if (any(!is.finite(x)))
    stop("missing values are not allowed")
  n <- length(x)
  if (n < 1L)
    stop("invalid length of 'x'") 
  xn <- (1:n)/n 
  lm.fit <- lm(x ~ 1 + xn + I(xn^2) + I(xn^3))
  denominator <- 4*(lm.fit[[1]][3])^2 + 12*lm.fit[[1]][3]*
    lm.fit[[1]][4] + 12*(lm.fit[[1]][4])^2
  q.n <- as.integer(floor((n)^(4/5)*(9/2)^(1/5)*
                            ((var(resid(lm.fit))/denominator)^(1/5))))
  qn <- ifelse(q.n < n, as.integer(min(q.n, n - q.n)), floor(n^(4/5)/2))
  return(qn)
}
