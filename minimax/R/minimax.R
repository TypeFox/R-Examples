pminimax <- function(x, a = 1, b = 1) {
  if(!is.numeric(x) | !is.numeric(a) | !is.numeric(b))
      stop("non-numeric arguments")
  ax <- attributes(x)
  n <- max(length(x), length(a), length(b))
  r <- rep(NA, n)
  x <- rep(x, length.out = n)
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)
  OK <- (a > 0 & b > 0)
  if(any(OK)) {
    x <- pmax(0, pmin(1, x[OK]))
    a <- a[OK]
    b <- b[OK]
    r[OK] <- 1 - (1-x^a)^b
  }
  attributes(r) <- ax
  r
}

dminimax <- function(x, a = 1, b = 1, log = FALSE) {
  if(!is.numeric(x) | !is.numeric(a) | !is.numeric(b))
      stop("non-numeric arguments")
  ax <- attributes(x)
  n <- max(length(x), length(a), length(b))
  r <- rep(NA, n)
  x <- rep(x, length.out = n)
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)
  OK <- (a > 0 & b > 0)
  if(any(OK)) {
    x <- pmax(0, pmin(1, x[OK]))
    a <- a[OK]
    b <- b[OK]
    if(log) 
        r[OK] <- log(a) + log(b) + (a-1)*log(x) + (b-1)*log(1 - x^a)
    else 
        r[OK] <- a*b*x^(a-1)*(1-x^a)^(b-1)
  }
  attributes(r) <- ax
  r
}

qminimax <- function(y, a = 1, b = 1) {
  if(!is.numeric(y) | !is.numeric(a) | !is.numeric(b))
      stop("non-numeric arguments")
  ay <- attributes(y)
  n <- max(length(y), length(a), length(b))
  r <- rep(NA, n)
  y <- rep(y, length.out = n)
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)
  OK <- (a > 0 & b > 0 & y >=0 & y <= 1)
  if(any(OK)) {
    y <- y[OK]
    a <- a[OK]
    b <- b[OK]
    r[OK] <- (1 - (1-y)^(1/b))^(1/a)
  }
  attributes(r) <- ay
  r
}

rminimax <- function(n, a = 1, b = 1) {
  an <- attributes(n)
  if(length(n) > 1)
      n <- length(n)
  if(!is.numeric(n) | !is.numeric(a) | !is.numeric(b))
     stop("non-numeric arguments")
  if(n < 0) stop("negative value for n")
  if(abs(n %% 1) > sqrt(.Machine$double.eps)) {
    warning("non-integer value for n")
    n <- round(n)
  }
  r <- rep(NA, n)
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)
  OK <- (a > 0 & b > 0)
  if(any(OK)) {
    a <- a[OK]
    b <- b[OK]
    n <- sum(OK)
    r[OK] <- (1 - runif(n)^(1/b))^(1/a)
  }
  attributes(r) <- an
  r
}
    
