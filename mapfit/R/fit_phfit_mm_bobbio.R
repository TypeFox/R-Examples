## Some methods from BuTools

hypoerlang <- function(shape,
  initprob = rep(1/length(shape),length(shape)),
  rate = rep(1,length(shape))) {
  size <- length(shape)
  phsize <- sum(shape)
  index <- cumsum(shape)
  sindex <- c(1, index + 1)[1:size]
  eindex <- index
  alpha <- numeric(phsize)
  xi <- numeric(phsize)
  alpha[sindex] <- initprob
  xi[phsize] <- rate[size]

  v <- numeric(0)
  i <- numeric(0)
  j <- numeric(0)
  for (k in 1:size) {
    i <- c(i, sindex[k]:eindex[k])
    j <- c(j, sindex[k]:eindex[k])
    v <- c(v, rep(-rate[k], shape[k]))
  }
  for (k in 1:size) {
    if (k != size) {
      i <- c(i, sindex[k]:eindex[k])
      j <- c(j, (sindex[k]+1):(eindex[k]+1))
      v <- c(v, rep(rate[k], shape[k]))
    } else {
      if (shape[k] != 1) {
        i <- c(i, sindex[k]:(eindex[k]-1))
        j <- c(j, (sindex[k]+1):eindex[k])
        v <- c(v, rep(rate[k], shape[k]-1))
      }
    }
  }
  Q <- sparseMatrix(dims=c(phsize,phsize), i=i, j=j, x=v)
  ph(alpha=alpha, Q=Q, xi=xi)
}

lowerbound <- function(n, n2) {
  if (n2<(n+1)/n) {
    lb <- Inf
  } else if (n2<(n+4)/(n+1)) {
    p <- ((n+1)*(n2-2)) / (3*n2*(n-1)) * ((-2*sqrt(n+1)) / sqrt(-3*n*n2+4*n+4) -1)
    a <- (n2-2) / (p*(1-n2) + sqrt(p^2+p*n*(n2-2)/(n-1)))
    l <- ((3+a)*(n-1)+2*a) / ((n-1)*(1+a*p)) - (2*a*(n+1)) / (2*(n-1)+a*p*(n*a+2*n-2))
    lb <- l
  } else {
    lb <- (n+1)/n *n2
  }
  lb
}

upperbound <- function(n, n2) {
  if (n2<(n+1)/n) {
    ub <- -Inf
  } else if (n2<=n/(n-1)) {
    u <- (2*(n-2)*(n*n2-n-1)*sqrt(1+(n*(n2-2))/(n-1)) + (n+2)*(3*n*n2-2*n-2)) / (n^2*n2)
    ub <- u
  } else {
    ub <- Inf
  }
  ub
}

isreal <- function(v) {
  if (abs(Im(v)) < sqrt(.Machine$double.eps)) {
    TRUE
  } else {
    FALSE
  }
}

mm.bobbio05 <- function(m1, m2, m3) {
  ## normalized moments
  n1 <- m1
  n2 <- m2 / m1^2
  n3 <- m3 / m1 / m2

  ## detect the number phases
  n <- 2
  while (n < 100 && ((n+1) / n > n2 || lowerbound(n, n2) >= n3 || upperbound(n, n2) <= n3)) {
    n <- n + 1
  }
  if (n2 < (n+1)/n) {
    n2 <- (n+1)/n
  }
  if (n3 < lowerbound(n, n2)) {
    n3 <- lowerbound(n, n2)
  }
  if (n3 > upperbound(n, n2)) {
    n3 <- upperbound(n, n2)
  }

  if (n2 > 2 || n3 < 2*n2 - 1) {
    b <- 2*(4-n*(3*n2-4)) / (n2*(4+n-n*n3) +
      sqrt(n*n2)*sqrt(12*n2^2*(n+1)+16*n3*(n+1)+n2*(n*(n3-15)*(n3+1)-8*(n3+3))))
    a <- (b*n2-2)*(n-1)*b / (b-1) / n
    p <- (b-1) / a
    lambda <- (p*a+1) / n1
    mu <- (n-1)*lambda / a
    res <- hypoerlang(shape=c(n-1,1), initprob=c(p, 1-p), rate=c(mu, lambda))
  } else {
    c4 <- n2*(3*n2-2*n3)*(n-1)^2;
    c3 <- 2*n2*(n3-3)*(n-1)^2;
    c2 <- 6*(n-1)*(n-n2);
    c1 <- 4*n*(2-n);
    c0 <- n*(n-2);
    fs <- polyroot(c(c0, c1, c2, c3, c4))
    found <- 0
    for (i in 1:length(fs)) {
      f <- fs[i]
      a <- 2*(f-1)*(n-1) / ((n-1)*(n2*f^2-2*f+2)-n)
      p <- (f-1)*a
      lambda <- (a+p) / n1
      mu <- (n-1) / (n1 - p/lambda)
      if (isreal(p) && isreal(lambda) && isreal(mu) && Re(p)>=0 && Re(p)<=1 && Re(lambda)>0 && Re(mu)>0) {
        p <- Re(p)
        lambda <- Re(lambda)
        mu <- Re(mu)
        res <- hypoerlang(shape=c(1,n-1), initprob=c(p, 1-p), rate=c(lambda, mu))
        found <- 1
      }
    }
    if (found == 0) {
      stop(sprintf("Cannot find the APH with 3 moments %f %f %f", m1, m2, m3))
    }
  }
  res
}
