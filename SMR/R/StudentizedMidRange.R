# computes the CDF of the  normal midrange
# input: the quantile q in R and  the number of means size >= 2. 
pNMR <- function(q, size, np = 32)
{
    if (q == -Inf) return(0.0) else if (q == Inf) return(1.0)    
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # interval 1: (-infty; q - 8]
    y  <- (q - 8) + (1 + x) / (x - 1) # from -1, 1 to -infty, q - 8
    aux <- pnorm(2*q-y)-pnorm(y)
    aux[aux <= 0] <- .Machine$double.eps^19 
    fy <- log(dnorm(y)) + (size-1) * log(aux)
    fy <- exp(fy)
    fy <- size * (2 / (x - 1)^2) * fy 
    I  <- sum(fy * w) 
    # interval 2: [q-8; q]
    a  <- q - 8
    b  <- q
    y  <- (b - a) / 2 * x + (a + b) / 2 # from -1, 1 to (q - 8),  q
    aux <- pnorm(2*q-y)-pnorm(y)
    aux[aux <= 0] <- .Machine$double.eps^19
    fy <- log(dnorm(y)) + (size-1) * log(aux)
    fy <- exp(fy)
    fy <- (b - a) / 2 * size * fy 
    I  <- I + sum(fy * w) 
    return(I) 
}

# auxiliar function 2
pMR_aux2 <- function(s, q, size, df, np = 32)
{
   
    Ii <- pNMR(s * q, size, np)    
    fx <- df/2 * log(df) - lgamma(df/2) - (df/2-1) * log(2) + 
                    (df - 1) * log(s) - df * s * s / 2
    fx <- exp(fx)
    I  <- fx * Ii
    return(I)  
}


# computes the CDF of the  externally studentized normal midrange
# input: the quantile q in R and  the number of means size >= 2,
# degrees of freedom df > 0, and number of points of the 
# gaussian quadrature (np)
pMR <- function(q, size, df, np = 32)
{
    if (is.nan(q) == TRUE) return(NaN)
    if (is.nan(size) == TRUE) return(NaN)
    if (is.nan(df) == TRUE) return(NaN)
    if (size<=1) return(NaN)
    if (df <= 0) return(NaN) 
    if (q == Inf) return(1)
    if (q == -Inf) return(0)
    if (df == Inf) return(pNMR(q, size, np))
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # Integration interval from 0 to 1
    y  <- -0.5 * x + 0.5 # from [-1; 1] to [0; 1]
    y  <- matrix(y, np, 1)
    fy <-  0.5 * apply(y, 1, pMR_aux2, q, size, df, np) 
    I  <- sum(fy * w)
    # integration interval from 1 to +infty
    y <- 1+(1+x)/(1-x) # from [-1; 1] to [1; +infty)
    y <- matrix(y, np, 1)
    fy.aux <- apply(y, 1, pMR_aux2, q, size, df, np) 
    fy <- log(fy.aux)+log(2)-2*log(1-x)
    fy <- exp(fy)
    I <- I+sum(fy * w) 
    return(I)
}

# computes the PDF of the normal midrange
# input: the quantile q in R and  the number of means size >= 2,
# degrees of freedom df > 0 and the number of points of the 
# gaussian quadrature (np). We use two divisions of the 
# integration limits. The first division was between Inf and q-8,
# and the second division was between q-8 and q.
dNMR <- function(q, size, np = 32)
{
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # interval 1: (-infty; q - 8]
    y  <- (q - 8) + (1 + x) / (x - 1) # from [-1; 1] to (-infty; q - 8]
    aux <- pnorm(2 * q - y) - pnorm(y)
    aux[aux <= 0] <- .Machine$double.eps^19
    fy <- log(dnorm(y)) + (size - 2) * log(aux) + log(dnorm(2 * q - y)) + log(2) +
          log(size) + log(size-1)
    fy <- exp(fy)
    fy <- (2 / (x - 1)^2) * fy 
    I  <- sum(fy * w) 
    # interval 2: [q-8; q]
    a  <- q - 8
    b  <- q
    y  <- (b - a) / 2 * x + (a + b) / 2 #from [-1; 1] to [(q - 8);  q]
    aux <- pnorm(2 * q - y) - pnorm(y)
    aux[aux <= 0] <- .Machine$double.eps^19 
    fy <- log(dnorm(y)) + (size - 2) * log(aux) + log(dnorm(2 * q - y)) + log(2) +
          log(size) + log(size-1)
    fy <- exp(fy)
    fy <- (b - a) / 2 * fy 
    I  <- I + sum(fy * w) 
    return(I) 
}

# auxiliar function 3
pMR_aux3 <- function(s, q, size, df, np = 32)
{   
    Ii <- dNMR(s * q, size, np)    
    fx <- df/2 * log(df) - lgamma(df/2) - (df/2-1) * log(2) + 
          (df - 1) * log(s) - df * s * s / 2
    fx <- exp(fx) * s
    I  <- fx * Ii
    return(I)  
}

# computes the PDF of the  externally studentized normal midrange
# input: the quantile q in R and the number of means size >= 2,
# degrees of freedom df > 0 and the number of points 
# of the gaussian quadrature (np). 
# We use two divisions of the integration limits. 
# The first division was between Inf and q-8, 
# and the second division was between q-8 and q.
dMR <- function(q, size, df, np = 32)
{
    if (is.nan(q) == TRUE) return(NaN)
    if (is.nan(size) == TRUE) return(NaN)
    if (is.nan(df) == TRUE) return(NaN)
    if (size <= 1) return(NaN)
    if (df <= 0) return(NaN) 
    if ((q == -Inf) | (q == Inf)) return(0)
    if (df == Inf) return(dNMR(q, size, np))
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # integration limit between 0 and 1
    y  <- -0.5 * x + 0.5 # changing the interval [-1,1] to [0,1]
    y  <- matrix(y, np, 1)
    fy <-  0.5 * apply(y, 1, pMR_aux3, q, size, df, np) 
    I  <- sum(fy * w) 
    # integration limit [1,Infty] - transformation y <- -ln(x), 0 < x < exp(-1)
    b <- exp(-1)
    a <- 0.0
    x <- (b - a) / 2.0 * x + (b + a) / 2.0 # changing the interval [-1,1] to [0,exp(0,-1)]        
    y  <- -log(x) # changing the interval [0,exp(-1)] to [1,infty]
    y  <- matrix(y, np, 1)
    fy <- (b - a) / 2.0 * apply(y, 1, pMR_aux3, q, size, df, np) / x 
    I  <- I + sum(fy * w) 
    return(I)  
}


# computes the quantiles of the normal midrange
# inputs: the cumulative probability (0 < p < 1), 
# the number of means (size >= 2) n and the number 
# of points of the gaussian quadrature (np)
qNMR <- function(p, size, np = 32, eps = 1.0e-13, maxit = 5000)
{
   if (p < 0.5) q0 <- -0.5 else
   if (p > 0.5) q0 <-  0.5 else q0 <- 0 # arbitrary initial value
   found <- FALSE   
   it <- 0   
   while ((found == FALSE) & (it <= maxit))
   {
      q1 <- q0 - (pNMR(q0, size, np) - p) / dNMR(q0, size, np)
      if (abs(q1-q0) <= eps) found = TRUE
      it <- it + 1
      q0 <- q1
   }
   return(q1)
}

# computes the quantiles of the  externally studentized normal midrange
# inputs: the cumulative probability (0 < p < 1), 
# the number of means (size >= 2) size, degrees of freedom df > 0 
# and number of points of the gaussian quadrature (np)
qMR <- function(p, size, df, np = 32, eps = 1e-13, maxit = 5000)
{
   if (is.nan(p) == TRUE) return(NaN)
   if (is.nan(size) == TRUE) return(NaN)
   if (is.nan(df) == TRUE) return(NaN)
   if (size <= 1) return(NaN)
   if (df <= 0) return(NaN)
   if (p == 1) return(Inf)
   if (p == 0) return(-Inf)   
   if (p < 0) return(NaN) else if (p > 1) return(NaN)   
   if (df == Inf) return(qNMR(p, size, np, eps, maxit)) 
   if (p < 0.5) q0 <- -0.5 else
   if (p > 0.5) q0 <-  0.5 else q0 <- 0 # arbitrary initial value
   found <- FALSE   
   it <- 0   
   while ((found == FALSE) & (it <= maxit))
   {
      q1 <- q0 - (pMR(q0, size, df, np) - p) / dMR(q0, size, df, np)
      if (abs(q1-q0) <= eps) found = TRUE
      it <- it + 1
      q0 <- q1
   }
   return(q1)
}

# computes the vector of random numbers of the normal midrange or
# externally studentized normal midrange.
# Inputs -  n: size of the vector to be simulated, size: sample size
# df in the interval [0, Inf], is the degrees of freedom.
rMR <- function(n, size, df = Inf)
{
    if (is.nan(n) == TRUE) return(NaN)
    if (is.nan(size) == TRUE) return(NaN)
    if (is.nan(df) == TRUE) return(NaN)
    if (n < 1) return(NaN)
    if (size <= 1) return(NaN)
    if (df <= 0) return(NaN)
    if (df == Inf) X <- (matrix(rnorm(n * size), n, size)) else
    X <- (matrix(rnorm(n * size), n, size)) / (rchisq(n, df) / df)^0.5 
    midrange <- function(x) 
    {
       return((max(x)+ min(x)) / 2)
    }
    res <- apply(X, 1, midrange)
    return(res)
}

