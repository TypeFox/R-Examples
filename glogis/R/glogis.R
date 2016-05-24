#######################################
## Generalized logistic distribution ##
#######################################

## cumulative distribution function of generalized logistic distribution
pglogis <- function(q, location = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)
{
  x <- (q - location)/scale
  if(log.p) {
    if(lower.tail) {
      shape * plogis(x, log.p = TRUE)
    } else {
      log(1 - plogis(x)^shape)
    }
  } else {
    if(lower.tail) {
      exp(shape * plogis(x, log.p = TRUE))
    } else {
      1 - plogis(x)^shape
    }
  }
}

## associated quantile function
qglogis <- function(p, location = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)
{
  q <- if(log.p) {
    if(lower.tail) {
      qlogis(p/shape, log.p = TRUE)
    } else {
      qlogis((-expm1(p))^(1/shape))
    }
  } else {
    if(lower.tail) {
      qlogis(log(p)/shape, log.p = TRUE)
    } else {
      qlogis((1 - p)^(1/shape))
    }
  }
  location + q * scale
}


## associated probbility density function
dglogis <- function(x, location = 0, scale = 1, shape = 1, log = FALSE)
{
  x <- (x - location)/scale
  pdf <- log(shape) - log(scale) - x - (shape + 1) * log(1 + exp(-x))
  if(log) pdf else exp(pdf)
}

## associated score function, i.e., gradient of dglogis(..., log = TRUE)
sglogis <- function(x, location = 0, scale = 1, shape = 1)
{
  ## scale observations
  x <- (x - location)/scale

  # derivatives of location
  rval1 <- 1/scale - (shape + 1) * (1/scale * exp(-x))/(1 + exp(-x))

  # derivatives of scale
  rval2 <- -1/scale + x/scale - (shape + 1) * ((x/scale) * exp(-x))/(1 + exp(-x))

  # derivatives of shape
  rval3 <- 1/shape - log(1 + exp(-x))

  rvalue <- cbind(rval1, rval2, rval3)
  colnames(rvalue) <- c("location", "scale", "shape")
  return(rvalue)
}

## random numbers via dumb inversion technique
rglogis <- function(n, location = 0, scale = 1, shape = 1) {
  qglogis(runif(n), location = location, scale = scale, shape = shape)
}
