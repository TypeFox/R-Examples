pgev <- function(q, loc = 0, scale = 1, shape = 0,
                 lower.tail = TRUE){
  if (min(scale) <= 0)
    stop("invalid scale parameter")

  if (length(shape) != 1)
    stop("'shape' should be a scalar")

  q <- (q - loc) / scale

  if (shape == 0)
    p <- exp(-exp(-q))

  else
    p <- exp(-pmax(1 + shape * q, 0)^(-1/shape))

  if (!lower.tail)
    p <- 1 - p

  return(p)
}

rgev <- function (n, loc = 0, scale = 1, shape = 0){
  
    if (min(scale) < 0) 
      stop("invalid scale parameter")
    
    if (length(shape) != 1) 
      stop("'shape' should be a scalar")
    
    if (shape == 0) 
      return(loc - scale * log(rexp(n)))
    
    else
      return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
}

dgev <- function (x, loc = 0, scale = 1, shape = 0,
                  log = FALSE){
  
  if (min(scale) <= 0) 
    stop("invalid scale parameter")
  
  if (length(shape) != 1) 
    stop("'shape' should be a scalar")
  
  x <- (x - loc)/scale

  if (shape == 0) 
    dns <- -log(scale) - x - exp(-x)

  else {
    n <- length(x)
    x <- 1 + shape * x
    xpos <- x[x > 0 | is.na(x)]
    scale <- rep(scale, length.out = n)[x > 0 | is.na(x)]
    dns <- rep(-Inf, n)
    dns[x > 0 | is.na(x)] <- -log(scale) - xpos^(-1/shape) - 
      (1/shape + 1) * log(xpos)
  }

  if (!log) 
    dns <- exp(dns)

  return(dns)
}

qgev <- function (p, loc = 0, scale = 1, shape = 0,
                  lower.tail = TRUE){
  
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
      1) 
    stop("`p' must contain probabilities in (0,1)")
  
  if (min(scale) < 0) 
    stop("invalid scale parameter")
  
  if (length(shape) != 1) 
    stop("'shape' should be a scalar")

  if (!lower.tail) 
    p <- 1 - p

  if (shape == 0) 
    return(loc - scale * log(-log(p)))
  
  else
    return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
}

rgpd <- function(n, loc = 0, scale = 1, shape = 0){
  
  if (min(scale) < 0) 
    stop("invalid scale parameter")
  
  if (length(shape) != 1) 
    stop("'shape' should be a scalar")
  
  if (shape == 0) 
    return(loc + scale * rexp(n))
  
  else
    return(loc + scale * (runif(n)^(-shape) - 1)/shape)
}

qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  lambda = 0){
  
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
      1) 
    stop("`p' must contain probabilities in (0,1)")
  
  if (min(scale) < 0) 
    stop("invalid scale parameter")
  
  if (length(shape) != 1) 
    stop("'shape' should be a scalar")
  
  if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1)
    stop("invalid lambda parameter")
  
  if (any(p < lambda))
    stop("``p'' must satisfy ``p >= lambda''")
  
  if (lower.tail) 
    p <- 1 - p
  
  p <- p / (1 - lambda)
  
  if (shape == 0) 
    return(loc - scale * log(p))
  
  else
    return(loc + scale * (p^(-shape) - 1)/shape)
}

dgpd <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE){
  
  if (min(scale) <= 0) 
    stop("invalid scale paramter")
  
  if (length(shape) != 1) 
    stop("'shape' should be a scalar")
  
  d <- (x - loc)/scale
  nn <- length(d)
  scale <- rep(scale, length.out = nn)
  index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)

  if (shape == 0) {
    d[index] <- -log(scale[index]) - d[index]
    d[!index] <- -Inf
  }

  else {
    d[index] <- -log(scale[index]) - (1/shape + 1) *
      log(1 + shape * d[index])
    d[!index] <- -Inf
  }
  
  if (!log) 
    d <- exp(d)

  return(d)
}

pgpd <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  lambda = 0){
  
  if (min(scale) <= 0) 
    stop("invalid scale parameter")
  
  if (length(shape) != 1) 
    stop("'shape' should be a scalar")
  
  if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1)
    stop("invalid lambda parameter")
  
  q <- pmax(q - loc, 0)/scale
  
  if (shape == 0) 
    p <- 1 - (1 - lambda) * exp(-q)
  
  else {
    p <- pmax(1 + shape * q, 0)
    p <- 1 - (1 - lambda) * p^(-1/shape)
  }
  
  if (!lower.tail) 
    p <- 1 - p
  
  return(p)
}
