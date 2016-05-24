
# Based on function "tseries::jarque.bera.test"

JarqueBera.test <- function(x, fc = 3.5, ...) UseMethod("JarqueBera.test")

JarqueBera.test.Arima <- function(x, fc = 3.5, ...)
{
  resid <- residuals(x)
  s <- frequency(x$residuals)

  id0 <- max(x$arma[c(1,6,2)]) + s * max(x$arma[c(3,7,4)])
  #id0 <- s * max(x$arma[c(3,7,4)])
  id0 <- which(resid[seq_len(id0)] == 0)
  if (length(id0) > 0)
    resid <- resid[-seq_len(max(id0))]
  tmp <- x$arma[6] + x$arma[5] * x$arma[7]
  id0 <- if (tmp > 1) seq.int(tmp) else c(1, 2)
  if (any(abs(resid[id0]) > fc * sd(resid[-id0])))
    resid <- resid[-id0]

  JarqueBera.test(resid)
}

JarqueBera.test.default <- function(x, ...)
{
  if (NCOL(x) > 1) 
    stop("x is not a vector or univariate time series")

  # do this before na.omit(x)
  DNAME <- deparse(substitute(x))

  #if (anyNA(x))) 
  #  stop("NAs in x")
  x <- na.omit(as.vector(x))

  n <- length(x)
  m1 <- sum(x)/n
  m2 <- sum((x - m1)^2)/n
  m3 <- sum((x - m1)^3)/n
  m4 <- sum((x - m1)^4)/n
  b1 <- (m3/m2^(3/2))^2
  b2 <- (m4/m2^2)
  STATISTIC <- n * b1/6 + n * (b2 - 3)^2/24  
  names(STATISTIC) <- "X-squared"
  #PVAL <- pchisq(STATISTIC, df = 2, lower.tail = FALSE)
  PVAL <- 1 - pchisq(STATISTIC, df = 2)

  sk <- sqrt(b1)
  skewness <- c(statistic = sk, 
    p.value = 2 * pnorm(abs(sk), mean = 0, sd = sqrt(6/n), lower.tail = FALSE))

  kurtosis <- c(statistic = b2, 
    p.value = 2 * pnorm(abs(b2 - 3), mean = 0, sd = sqrt(24/n), lower.tail = FALSE))

  jb <- structure(list(statistic = STATISTIC, parameter = c(df = 2), 
    p.value = PVAL, method = "Jarque Bera Test", data.name = DNAME), 
    class = "htest")   

  skewness <- structure(list(statistic = skewness[1], parameter = NULL, 
    p.value = skewness[2], method = "Skewness", data.name = DNAME), 
    class = "htest")
  
  kurtosis <- structure(list(statistic = kurtosis[1], parameter = NULL, 
    p.value = kurtosis[2], method = "Kurtosis", data.name = DNAME), 
    class = "htest")
    
  structure(list(jb, skewness, kurtosis), class = "mhtest")
}

print.mhtest <- function(x, digits = 4L, ...)
{
  # "x" is a list containing "htest" objects
  for (i in seq_along(x))
  {
    stopifnot(inherits(x[[i]], "htest"))
    print(x[[i]], ...)   
  }
}
