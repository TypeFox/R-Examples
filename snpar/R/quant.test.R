quant.test <-
function(x, y = NULL, q, paired = FALSE, p = 0.5, 
                      alternative = c("two.sided","less","greater"), exact = FALSE,
                      correct = TRUE)
{
  DNAME <- deparse(substitute(x))
  alternative <- match.arg(alternative)
  if (is.null(y))
  {
    METHOD <- ifelse(exact,"Exact one-sample quantile test",
                     "Approximate one-sample quantile test")
    if (!missing(q) && ((length(q) > 1L) || !is.finite(q))) 
      stop("'q' must be a finite single number")
    if (paired)
      stop("argument 'paired' can be only used for two-sample test")
  }  
  else
  {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    if (length(unique(y)) == 1)
      stop("please use argument 'q' for one-sample test")
    else if (length(unique(y)) == 2)
    {
      DNAME <- paste(DNAME, "with group", deparse(substitute(y)))
      xx <- x[y == unique(y)[1]]; yy <- x[y == unique(y)[2]]
      x <- NULL; y <- NULL; x <- xx; y <- yy
    }
    else
      DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    y <- y[is.finite(y)]
    n <- as.integer(length(y))
    if (n < 2L) 
      stop(paste("not enough 'y' data or please specify q =", y, "for one-sample test"))
    if (paired)
    {
      if (length(x) != n)
        stop("'x' and 'y' must have the same length")
      OK <- complete.cases(x,y)
      x <- x[OK] - y[OK]
      METHOD <- ifelse(exact,"Exactly paired two-sample quantile test",
                       "Approximately paired two-sample quantile test")
    }
    else
      METHOD <- ifelse(exact,"Exact two-sample (Brown-Mood) quantile test",
                       "Approximate two-sample (Brown-Mood) quantile test")
  }
  x <- x[is.finite(x)]
  m <- as.integer(length(x))
  if (m < 1L) 
    stop("not enough 'x' data")
  if (!is.numeric(x))
    stop("'x' must be numeric") 
  if (!(length(p) == 1L  && (p > 0) && (p < 1)))
    stop("'p' must be a single number between 0 and 1")  
  if (is.null(y) || (!is.null(y) && paired))
  { 
    q <- if (!is.null(y) && paired) 0 else q
    x <- x - q
    zeros <- any(x == 0)
    if (zeros) 
      x <- x[x != 0]
    s.p <- sum(x > 0)
    s.m <- m - s.p
    p.val.2sided <- ifelse(exact,2*min(pbinom(min(s.p,s.m),m,p), 1 - pbinom(min(s.p,s.m),
                           m,p)), 2*min(pnorm(min(s.p,s.m), m*p,sqrt(m*p*(1-p))),
                            1 - pnorm(min(s.p,s.m), m*p,sqrt(m*p*(1-p)))))  
    p.val.lsided <- ifelse(exact,pbinom(s.p,m,p),pnorm(s.p,m*p,sqrt(m*p*(1-p))))
    p.val.gsided <- ifelse(exact,pbinom(s.m-1,m,p),pnorm(s.m,m*p,sqrt(m*p*(1-p))))
    p.val <- switch(alternative, two.sided = p.val.2sided, less = p.val.lsided,
                    greater = p.val.gsided)
    ESTIMATE <- setNames(quantile(x,p)[[1]] + q,"location")
    alt <- switch(alternative,two.sided = paste("true location is not equal to", q),
                  less = paste("true location is less than", q), 
                  greater = paste("true location is greater than", q))
    result <- list(statistic = ESTIMATE, p.value = p.val, alternative = alt)
  }
  else
  {
    quanEST <- as.numeric(quantile(c(x,y),p))
    ESTIMATE <- setNames(as.numeric(quantile(x,p) - quantile(y,p)),"Difference")
    a <- sum(x > quanEST); b <- sum(y > quanEST)
    if (correct)
    {
      if (a < m/2)
        Z <- (a + 0.5 - m*(a + b)/(m + n))/sqrt(m*n*(a + b)*(m + n - a - b)/(m + n)^3)
      else
        Z <- (a - 0.5 - m*(a + b)/(m + n))/sqrt(m*n*(a + b)*(m + n - a - b)/(m + n)^3)
    }
    else
      Z <- (a - m*(a + b)/(m + n))/sqrt(m*n*(a + b)*(m + n - a - b)/(m + n)^3)
    p.val.2side <- ifelse(exact,2*min(phyper(a,m,n,a+b),1 - phyper(a,m,n,a+b)),
                          2*min(pnorm(Z), 1 - pnorm(Z)))
    p.val.lside <- ifelse(exact, phyper(a,m,n,a+b), pnorm(Z))
    p.val.gside <- ifelse(exact, 1 - phyper(a-1,m,n,a+b), 1 - pnorm(Z))
    p.val <- switch(alternative, two.sided = p.val.2side, less = p.val.lside,
                    greater = p.val.gside)
    xname <- deparse(substitute(x)); yname <- deparse(substitute(y))
    alt <- switch(alternative,two.sided = "true location shift is not equal to 0",
                  less = "true location shift is less than 0",
                  greater = "true location shift is greater than 0")
    result <- list(statistic = ESTIMATE, p.value = p.val, alternative = alt)
  }
  RVAL <- c(list(data.name = DNAME, method = METHOD), result)
  class(RVAL) <- "htest"
  return(RVAL)
}
