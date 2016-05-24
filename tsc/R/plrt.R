plrt <-
function(x, y, k = 500)
{
  # x is the vector containg the x sample values
  # y is the vector containg the y sample values
  # k is the number of Legendre quadrature points
  
  # The output of plrt includes the value lam of the likelihood ratio test statistic 
  # for testing the equality of two independent normal distributions, and the P-value
  # of the test P(lambda_{n,m} <= lam)
  
  # Generate the Legendre quadrature points
  # r contains the nodes, w contains the weights
  Legendre <- function(p) {
    x <- matrix(0, p, p)
    b <- c(1:(p - 1))^2
    b <- sqrt(b/(4 * b - 1))
    for(i in 2:p) {
      x[i, i - 1] <- b[i - 1]
    }
    x <- x + t(x)
    a <- eigen(x)
    u <- a$values
    v <- a$vectors
    w <- 2 * v[1,  ]^2
    cbind(u[p:1], w[p:1])
  }
  z <- Legendre(k)
  r <- z[, 1]
  w <- z[, 2]
  
  # Compute the likelihood ratio test statistic
  n <- length(x)
  m <- length(y)
  xb <- mean(x)
  yb <- mean(y)
  u <- mean(c(x, y))
  num <- log(sum((x - xb)^2)/n)*(n/2) +log(sum((y - yb)^2)/m)*(m/2)
  den <- log((sum((x - u)^2) + sum((y - u)^2))/(n + m))*((n + m)/2)
  lam1<- num-den
  lam<-exp(lam1)
  
  
  # Locate the intervals containing the two roots a and b, and find a and b
  u <- seq(0,1,0.005)
  v <- (exp(log(lam)+n/2*log(n)+m/2*log(m)-(n+m)/2*log(n+m)-n/2*log(u)))^(2/m)
  u <- u[1-u>v]
  fn <- function(w, ...)
  {
    1-w - (exp(log(lam)+n/2*log(n)+m/2*log(m)-(n+m)/2*log(n+m)-n/2*log(w)))^(2/m)
  }
  
  fn2 <- function(w, ...) fn(w)^2
  lower = 1e-10
  upper = u[1]
  if (fn(lower)*fn(upper)<=0) a <- uniroot(fn, lower = 1e-10, upper = u[1], tol=1e-8, lam = lam, n = n, m = m)$root
  if (fn(lower)*fn(upper)>0) a<-optimize(fn2,c(lower,upper),tol=0.00001)$minimum
  lower = u[1]
  upper = 1-1e-10
  if (fn(lower)*fn(upper)<=0) b <- uniroot(fn, lower = u[1], upper = 1-1e-10, tol=1e-8, lam = lam, n = n, m = m)$root
  if (fn(lower)*fn(upper)>0) b<-optimize(fn2,c(lower,upper),tol=0.00001)$minimum
  # Compute the double integral
  h1 <- (b - a)/2
  h2 <- (b + a)/2
  w1 <- h1 * r + h2
  d1 <- 1 - w1
  c1 <- (exp(log(lam)+n/2*log(n)+m/2*log(m)-(n+m)/2*log(n+m)-n/2*log(w1)))^(2/m)
  k1 <- (d1 - c1)/2
  k2 <- (d1 + c1)/2
  s <- rep(0,k)
  for(i in 1:k) {
    w2 <- k1[i] * r + k2[i]
    s[i] <- k1[i] * sum((w1[i]^((n - 1)/2 - 1) * w2^((m - 1)/2 - 1))/sqrt(1 - w1[i] - w2) * w)
  }                            
  v <- h1 * sum(s*w)
  pv <- 1 - exp(lgamma((n + m - 1)/2) - lgamma((n - 1)/2) - lgamma((m - 1)/2) - lgamma(0.5))*v
  final <- list(test_stat=lam, p_value=pv)
  return(final)
  
  # cat("Likelihood Ratio Test Statistic = ", round(lam, 4), "\n")
  #cat("P-value = ", round(pv, 4), "\n")
}
