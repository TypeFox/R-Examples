# generates random sample from multivariate normal distribution
# with means mu and covariance sigma positive definite
# input: the sample size N, the parameters mu and sigma
rmultvariate <- function (N = 1, mu, Sigma, tol = 1e-06) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")  
  Sd <- eigen(Sigma, symmetric = TRUE)
  eva <- Sd$values
  if (!all(eva >= -tol * abs(eva[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * N), N)  
  X <- drop(mu) + Sd$vectors %*% diag(sqrt(pmax(eva, 0)), p) %*% t(X)  
  if (N == 1) 
    drop(X) else t(X)
}


# Function to calculate nodes and weights of Gauss-Legendre 
# quadrature, where n is the number of points. It uses the  
# Golub-Welsh algorithm                                     
GaussLegendre <- function(n) 
{
  n <- as.integer(n)
  if (n < 0) 
    stop("Must be a non-negative number of nodes!")
  if (n == 0) 
    return(list(x = numeric(0), w = numeric(0)))
  i  <- 1:n
  j   <- 1:(n-1)
  mu0 <- 2
  b <- j / (4 * j^2 - 1)^0.5
  A <- rep(0, n * n)
  A[(n + 1) * (j - 1) + 2] <- b
  A[(n + 1) * j] <- b
  dim(A) <- c(n, n)
  sd <- eigen(A, symmetric = TRUE)
  w <- rev(as.vector(sd$vectors[1, ]))
  w <- mu0 * w^2
  x <- rev(sd$values)
  return(list(nodes = x, weights = w))
}

# Function to compute cumulative distribution  function of the noncentral 
# unilateral Dunnett's test statistic with infinity degrees of freedom nu
# given the vector of quantile q, the vector of correlation r, the vector 
# of the noncentrality parameter delta, the vector of degrees of  freedom
# nu, the number of points n of  the Gauss-Legendre  quadrature  and  the 
# list x of nodes and weights of the quadrature. It returns    cumulative 
# probabilities.
pNDUD <- function(q, r, delta, n = 32, x)
{
   k <- length(r)
   # transformando de -1 a 1 para -infty a +infty, t = x/(1-x^2), x em -1 a 1
   y <- matrix(x$nodes/(1-(x$nodes)^2), n, 1)   
   pArgYi <- function(y, r, delta)
   {
     ti <- pnorm((r^0.5 * y + q - delta)/(1-r)^0.5 , log.p=TRUE) 
     ti <- exp(sum(ti)+dnorm(y, log=TRUE))
     return(ti)
   }
   I <- apply(y, 1, pArgYi, r, delta)
   I <- I * ((1+x$nodes^2)/(1-(x$nodes)^2)^2)
   I <- sum(I * x$weights)
   return(I)
}

# Function to compute cumulative distribution  function of the noncentral 
# bilateral  Dunnett's test statistic with infinity degrees of freedom nu
# given the vector of quantile q, the vector of correlation r, the vector 
# of the noncentrality parameter delta, the vector of degrees of  freedom
# nu, the number of points n of  the Gauss-Legendre  quadrature  and  the 
# list x of nodes and weights of the quadrature. It returns    cumulative 
# probabilities.
pNDBD <- function(q, r, delta, n = 32, x)
{
   k <- length(r)
   # transform from [-1; 1] to (-infty, +infty), t = x/(1-x^2), x in [-1, 1]
   y <- matrix(x$nodes/(1-(x$nodes)^2), n, 1)   
   pArgYi <- function(y, r, delta)
   {
     ti <- pnorm((r^0.5 * y + q - delta)/(1-r)^0.5 , log.p=FALSE) 
     ti <- ti - pnorm((r^0.5 * y - q - delta)/(1-r)^0.5, log.p=FALSE)
     if (any(ti <= 0)) 
     {
        ti[ti <= 0] <- -743.7469 
        ti[ti > 0] <- log(ti[ti > 0])
     } else ti <- log(ti) 
     ti <- exp(sum(ti)+dnorm(y, log=TRUE))
     return(ti)
   }
   I <- apply(y, 1, pArgYi, r, delta)
   I <- I * ((1+x$nodes^2)/(1-x$nodes^2)^2)
   I <- sum(I * x$weights)
   return(I)
}


# Function to compute the probability density function of the  noncentral 
# unilateral Dunnett's test statistic with infinity degrees of freedom nu
# given the vector of quantile q, the vector of correlation r, the vector 
# of the noncentrality parameter delta, the vector of degrees of  freedom
# nu, the number of points n of  the Gauss-Legendre  quadrature  and  the 
# list x of nodes and weights of the quadrature.   It  returns    density 
# values.
dNDUD <- function(q, r, delta, n = 32, x)
{
  k <- length(r)
  # transform from [-1; 1] to (-infty, +infty), t = x/(1-x^2), x in [-1, 1]
  y <- matrix(x$nodes/(1-(x$nodes)^2), n, 1)  
  rgdgi <- cbind(r, delta, (1:k))
  loopg <-  function(rd, ym, q)
  {
    tg <- pnorm((rd[1]^0.5 * ym + q - rd[2]) / (1 - rd[1])^0.5, log.p=TRUE) 
    return(tg) 
  }    
  loopi <- function(rd, ym, q)
  {
    ti <- dnorm((rd[1]^0.5 * ym + q - rd[2]) / (1 - rd[1])^0.5, log=TRUE) - 0.5 * log(1-rd[1])     
    return(ti)
  }    
  # loop of y's
  loopy <- function(ym, q, rd)
  {
    Phy <- apply(rd, 1, loopg, ym, q) # loopg
    phy <- apply(rd, 1, loopi, ym, q) # loopi
    sPhy <- sum(Phy) # sum of  Phi's differences
    sTi <- sum(exp(phy + sPhy - Phy)) * dnorm(ym)
    return(sTi)
  }  
  I <- apply(y, 1, loopy, q, rgdgi)
  I <- I * ((1+x$nodes^2)/(1-x$nodes^2)^2)
  I <- sum(I * x$weights) 
  return(I)
}

# Function to compute the probability density function of the  noncentral 
# bilateral Dunnett's test statistic with infinity degrees of freedom  nu
# given the vector of quantile q, the vector of correlation r,the vector 
# of the noncentrality parameter delta, the vector of degrees  of freedom
# nu, the number of points n of  the Gauss-Legendre  quadrature  and  the 
# list x of nodes and  weights of the  quadrature.    It  returns density
# probabilities.
dNDBD <- function(q, r, delta, n = 32, x)
{
  k <- length(r)
  # transform from [-1; 1] to (-infty, +infty), t = x/(1-x^2), x in [-1, 1]
  y <- matrix(x$nodes/(1-(x$nodes)^2), n, 1)  
  rgdgi <- cbind(r, delta, (1:k))
  loopg <-  function(rd, ym, q)
  {
      tg <- pnorm((rd[1]^0.5 * ym + q - rd[2]) / (1 - rd[1])^0.5) - 
            pnorm((rd[1]^0.5 * ym - q - rd[2]) / (1 - rd[1])^0.5)
      if (tg <= 0) tg <- -743.7469 else tg <- log(tg) 
      return(tg) 
  }    
  loopi <- function(rd, ym, q)
  {
      ti <- dnorm((rd[1]^0.5 * ym + q - rd[2]) / (1 - rd[1])^0.5) + 
            dnorm((rd[1]^0.5 * ym - q - rd[2]) / (1 - rd[1])^0.5)
      if (ti <= 0) ti <- -743.7469 else ti <- log(ti) 
      aux <-  ti - 0.5 * log(1 - rd[1])
      return(aux)
  }    
  # loop of y's
  loopy <- function(ym, q, rd)
  {
     Phy <- apply(rd, 1, loopg, ym, q) # loopg
     phy <- apply(rd, 1, loopi, ym, q) # loopi
     sPhy <- sum(Phy) # sum of Phi's differences
     sTi <- sum(exp(phy + sPhy - Phy)) * dnorm(ym)
     return(sTi)
  }  
  I <- apply(y, 1, loopy, q, rgdgi)
  I <- I * ((1+x$nodes^2)/(1-x$nodes^2)^2)
  I <- sum(I * x$weights) 
  return(I)
}

# Function to compute      quantiles   from   the  noncentral unilateral 
# Dunnett's test statistic distribution with infinity degrees of freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and weights  of  the  quadrature. It  returns   
# quantiles.
qNDUD <- function(p, r, delta, n = 32, x)
{
   k <- length(r)
   alpha <- 1 - p
   alpha_k <- 1 - (1-alpha)^(1/k)
   q0 <- qnorm(1 - alpha_k) # initial value
   maxIt <- 5000
   tol <- 1e-13
   conv <- FALSE
   it <- 0
   while ((conv == FALSE) & (it <= maxIt))
   {
       q1 <- q0 - (pNDUD(q0, r, delta, n, x) - p) / dNDUD(q0, r, delta, n, x)
       if (abs(q1-q0) <= tol) conv <- TRUE
       q0 <- q1
       it <- it + 1
   }
   return(q1)  
}

# Function to compute      quantiles  from    the   noncentral bilateral 
# Dunnett's test statistic distribution with infinity degrees of freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and  weights of the quadrature.  It   returns 
# quantiles.
qNDBD <- function(p, r, delta, n = 32, x)
{
   k <- length(r)
   alpha <- 1 - p
   alpha_k <- 1 - (1-alpha/2)^(1/k)
   q0 <- qnorm(1 - alpha_k) # initial value
   if (q0 < 0) q0 <- -q0
   maxIt <- 5000
   tol <- 1e-13
   conv <- FALSE
   it <- 0
   while ((conv == FALSE) & (it <= maxIt))
   {
       q1 <- q0 - (pNDBD(q0, r, delta, n, x) - p) / dNDBD(q0, r, delta, n, x)
       if (abs(q1-q0) <= tol) conv <- TRUE
       q0 <- q1
       it <- it + 1
   }
   return(q1)  
}

# Function to compute the Cumulative  distribution   of  the  noncentral  
# unilateral Dunnett's test statistic with finite   degrees  of  freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and weights  of  the quadrature.   It returns   
# cumulative probabilities.
pNDUDF <- function(q, r, nu, delta, n = 32, xx)
{
  k <- length(r)
  x  <- xx$nodes
  w  <- xx$weights
  # computing integral in [0,1]
  y  <- 0.5 * x + 0.5 # from [-1, 1] to [0, 1]
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + (nu-1) * log(y) - nu * y * y / 2
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <-  0.5 * apply(y, 1, pNDUD, r, delta, n, xx) * fx
  I  <- sum(fy * w) 
  # computing integral in [1, infty), using the transformation 
  # y <- -ln(x), 0 < x < exp(-1)  
  b <- exp(-1)
  a <- 0
  x <- (b - a) / 2 * x + (b + a) / 2 # from [-1,1] to [0, exp(-1)]
  y  <- -log(x) # from [0, exp(-1)] to [1, +infty)
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + (nu-1) * log(y) - nu * y * y / 2
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <- (b - a) / 2 * apply(y, 1, pNDUD, r, delta, n, xx) / x * fx
  I  <- I + sum(fy * w) 
  return(I)  
}

# Function to compute the Cumulative  distribution   of  the  noncentral  
# bilateral  Dunnett's test statistic with finite   degrees  of  freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and   weights of the quadrature.  It  returns   
# cumulative probabilities.
pNDBDF <- function(q, r, nu, delta, n = 32, xx)
{
  k <- length(r)
  x  <- xx$nodes
  w  <- xx$weights
  # computing integral in [0,1]
  y  <- 0.5 * x + 0.5 # from [-1,1] to [0, 1]
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + (nu-1) * log(y) - nu * y * y / 2
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <-  0.5 * apply(y, 1, pNDBD, r, delta, n, xx) * fx
  I  <- sum(fy * w) 
  # computing integral in [1, infty), using the transformation 
  # y <- -ln(x), 0 < x < exp(-1)  
  b <- exp(-1)
  a <- 0
  x <- (b - a) / 2 * x + (b + a) / 2 # from [-1,1] to [0, exp(-1)]        
  y  <- -log(x) # from [0, exp(-1)] to [1, +infty)
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + (nu-1) * log(y) - nu * y * y / 2
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <- (b - a) / 2 * apply(y, 1, pNDBD, r, delta, n, xx) / x * fx
  I  <- I + sum(fy * w) 
  return(I)  
}

# Function to compute the probability density function of the noncentral  
# unilateral Dunnett's test statistic with finite   degrees  of  freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and weights of the quadrature.   It   returns   
# density values.
dNDUDF <- function(q, r, nu, delta, n = 32, xx)
{
  k <- length(r)
  x  <- xx$nodes
  w  <- xx$weights
  # computing integral in [0,1]
  y  <- 0.5 * x + 0.5 # from [-1,1] to [0, 1]
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + 
            (nu-1) * log(y) - nu * y * y / 2 + log(y) # + jacobian
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <-  0.5 * apply(y, 1, dNDUD, r, delta, n, xx) * fx
  I  <- sum(fy * w) 
  # computing integral in [1, infty), using the transformation 
  # y <- -ln(x), 0 < x < exp(-1)  
  b <- exp(-1)
  a <- 0
  x <- (b - a) / 2 * x + (b + a) / 2 # from [-1,1] to [0, exp(-1)] 
  y  <- -log(x) # from [0, exp(-1)] to [1, +infty)
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + 
          (nu-1) * log(y) - nu * y * y / 2  + log(y) # + jacobian
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <- (b - a) / 2 * apply(y, 1, dNDUD, r, delta, n, xx) / x * fx
  I  <- I + sum(fy * w) 
  return(I)  
}

# Function to compute the probability density function of the noncentral  
# bilateral  Dunnett's test statistic with finite   degrees  of  freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and weights of the quadrature.   It   returns   
# density values.
dNDBDF <- function(q, r, nu, delta, n = 32, xx)
{
  k <- length(r)
  x  <- xx$nodes
  w  <- xx$weights
  # computing integral in [0,1]
  y  <- 0.5 * x + 0.5 # from [-1,1] to [0, 1]
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + 
           (nu-1) * log(y) - nu * y * y / 2 + log(y) # + jacobian
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <-  0.5 * apply(y, 1, dNDBD, r, delta, n, xx) * fx
  I  <- sum(fy * w) 
  # computing integral in [1, infty), using the transformation 
  # y <- -ln(x), 0 < x < exp(-1) 
  b <- exp(-1)
  a <- 0
  x <- (b - a) / 2 * x + (b + a) / 2 # from [-1,1] to [0, exp(-1)]        
  y  <- -log(x) # from [0, exp(-1)] to [1, +infty)
  fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + 
         (nu-1) * log(y) - nu * y * y / 2 + log(y) # + jacobian
  fx <- exp(fx)
  y  <- matrix(y, n, 1) * q
  fy <- (b - a) / 2 * apply(y, 1, dNDBD, r, delta, n, xx) / x * fx
  I  <- I + sum(fy * w) 
  return(I)  
}
 
# Function    to    compute      quantiles     from    the    noncentral  
# unilateral Dunnett's test statistic with finite   degrees  of  freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes  and weights  of the quadrature.   It returns  
# quantiles.
qNDUDF <- function(p, r, nu, delta, n = 32, x)
{
  k <- length(r)
  alpha <- 1 - p
  alpha_k <- 1 - (1-alpha)^(1/k)
  q0 <- qt(1 - alpha_k, nu) + mean(delta)# initial value 
  p0 <- pNDUDF(q0, r, nu, delta, n, x)
  if (abs(p0-p) > 0.1) q0  <- qNDUD(p, r, delta, n, x)
  maxIt <- 5000
  tol <- 1e-11
  conv <- FALSE
  it <- 0
  while ((conv == FALSE) & (it <= maxIt))
  {
    q1 <- q0 - (pNDUDF(q0, r, nu, delta, n, x) - p) / dNDUDF(q0, r, nu, delta, n, x)
    if (abs(q1-q0) <= tol) conv <- TRUE
    if (q1 < 0) q1 <- -q1
    q0 <- q1
    it <- it + 1
  }
  return(q1)  
}

# Function  to    compute         quantiles from    the       noncentral  
# bilateral  Dunnett's test statistic with finite   degrees  of  freedom 
# nu given the vector of quantile q, the vector of   correlation r,  the  
# vector of the noncentrality parameter delta, the vector of degrees  of 
# freedom nu, the number of points n of  the Gauss-Legendre   quadrature  
# and  the list x of nodes and weights of the quadrature.   It   returns   
# quantiles.
qNDBDF <- function(p, r, nu, delta, n = 32, x)
{
  k <- length(r)
  alpha <- 1 - p
  alphak <- 1 - (1-alpha/2)^(1/k)
  q0 <- qt(1 - alphak, nu) + mean(delta) # valor 
  if (q0 < 0) q0 <- -q0
  p0 <- pNDBDF(q0, r, nu, delta, n, x)
  if (abs(p0-p) > 0.1) q0  <- qNDBD(p, r, delta, n, x)
  maxIt <- 5000
  tol <- 1e-11
  conv <- FALSE
  it <- 0
  while ((conv == FALSE) & (it <= maxIt))
  {
    p0 <- pNDBDF(q0, r, nu, delta, n, x)
    q1 <- q0 - (p0 - p) / dNDBDF(q0, r, nu, delta, n, x)
    if (abs(q1-q0) <= abs(q0) * tol) conv <- TRUE
    if (q1 < 0) q1 <- -q1
    q0 <- q1
    it <- it + 1
  }
  return(q1)  
}