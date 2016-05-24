runs.test <-
function(x, exact = FALSE, alternative = c("two.sided", "less", "greater"))
{
  DNAME <- deparse(substitute(x))
  METHOD <- ifelse(exact,"Exact runs test", "Approximate runs rest")
  alternative <- match.arg(alternative)
  x <- x[is.finite(x)]
  N <- as.integer(length(x))
  if (N < 1L) 
    stop("not enough 'x' data")
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  x.s <- if (length(unique(x)) == 2) (x == unique(x)[1])*1 
  else (x > median(as.matrix(x)))*1
  m <- sum(1 - x.s)
  n <- N - m
  R <- 1
  for (i in 1:(N - 1))
  {
    if (x.s[i] != x.s[i+1])  R <- R + 1
  }
  STATISTIC <- setNames(R, "Runs")
  P.even <- function(m,n,r) 2*choose(m-1,r-1)*choose(n-1,r-1)/choose(m + n,n)
  P.odd <- function(m,n,r)
    (choose(m-1,r-1)*choose(n-1,r) + choose(n-1,r-1)*choose(m-1,r))/choose(m + n,n) 
  if (exact)
  {
    if (any(is.na(P.even(m,n,1:floor(R/2)))) || any(is.na(P.odd(m,n,1:floor(R/2)))))
      stop("can't calculate exact p-value; please use approximate method")
    if (R%%2 == 0)
    {
      p.val <- sum(P.even(m,n,1:(R/2))) + sum(P.odd(m,n,1:(R/2 - 1)))
      p.val1 <- 1 - p.val + P.even(m,n,R/2)
    }
    else 
    {
      p.val <- sum(P.even(m,n,1:floor(R/2))) + sum(P.odd(m,n,1:floor(R/2))) 
      p.val1 <- 1 - p.val + P.odd(m,n,floor(R/2))
    }        
  }
  else
    Z <- (R - 2*m*n/N - 1)/sqrt(2*m*n*(2*m*n - N)/(N^2*(N - 1)))
  P.VAL <- switch(alternative, two.sided = ifelse(exact, 2*min(p.val,1 - p.val), 
                  2*min(pnorm(Z),1- pnorm(Z))), less = ifelse(exact, p.val, pnorm(Z)),
                  greater = ifelse(exact, p.val1, 1 - pnorm(Z)))
  RVAL <- list(data.name = DNAME, method = METHOD, alternative = alternative,
               statistic = STATISTIC, p.value = P.VAL)
  class(RVAL) <- "htest"
  return(RVAL)
}
