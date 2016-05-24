#############################################################################
##
##  Transformed Density.
##
#############################################################################

Tfdd <- function(lpdf, dlpdf, d2lpdf, cT, x) {
  ## -----------------------------------------------------------------------
  ## Compute transformed density and its derivatives by means of
  ## log-density and its derivatives.
  ## -----------------------------------------------------------------------
  ##   lpdf   ... log-density 
  ##   dlpdf  ... derivative of log-density
  ##   d2lpdf ... 2nd derivative of log-density
  ##   cT     ... parameter for transformation
  ##   x      ... argument
  ## -----------------------------------------------------------------------
  ## Return: vector consisting of
  ##   Tfx    ... transformed density at x
  ##   dTfx   ... derivative of transformed density at x
  ##   d2Tfx  ... second derivative of transformed density at x
  ## -----------------------------------------------------------------------

  ## Remove attributes (s.t. 'identical' works as expected).
  cT <- as.numeric(cT)
  
  ## Evaluate log density and its derivatives.
  lfx <- lpdf(x)
  dlfx <- dlpdf(x)
  d2lfx <- d2lpdf(x)

  if (identical(cT, 0)) {
    ## Case: T(x) = log(x)
    return (c(lfx,dlfx,d2lfx))
  }

  else {
    ## Case: T(x) = sign(c) * x^c
    Tfx   <- sign(cT) * exp(cT * lfx)
    dTfx  <- cT * Tfx * dlfx
    d2Tfx <- cT * Tfx * (cT * dlfx^2 + d2lfx)
    return (c(Tfx,dTfx,d2Tfx))
  }
}

## --------------------------------------------------------------------------

Tf <- function(lpdf, cT, x) {
  ## -----------------------------------------------------------------------
  ## Compute transformed density by means of its log-density.
  ## -----------------------------------------------------------------------
  ##   lpdf   ... log-density 
  ##   cT     ... parameter for transformation
  ##   x      ... argument
  ## -----------------------------------------------------------------------
  ## Return: transformed density at x.
  ## -----------------------------------------------------------------------

  ## Remove attributes (s.t. 'identical' works as expected).
  cT <- as.numeric(cT)

  return (if (identical(cT, 0)) {lpdf(x)} else {sign(cT) * exp(cT * lpdf(x))})
}

## --------------------------------------------------------------------------

Tinv <- function(cT, x) {
  ## --------------------------------------------------------------------------
  ## Compute inverse transformation.
  ## -----------------------------------------------------------------------
  ##   cT     ... parameter for transformation
  ##   x      ... argument
  ## -----------------------------------------------------------------------
  ## Return: inverse transformation at x.
  ## -----------------------------------------------------------------------

  ## Remove attributes (s.t. 'identical' works as expected).
  cT <- as.numeric(cT)

  if (identical(cT, 0))
    ## Case: T(x) = log(x)
    return (exp(x))

  if (identical(cT, -0.5))
    ## Case: T(x) = -1/sqrt(x)
    return (1/(x*x))
  
  if (identical(cT, 1))
    ## Case: T(x) = x
    return (x)
  
  else
    ## Case: T(x) = sign(c) * x^c
    return ((sign(cT)*x)^(1/cT))
}

## --------------------------------------------------------------------------

FT <- function(cT, x) {
  ## --------------------------------------------------------------------------
  ## Compute antiderivative of inverse transformation.
  ## -----------------------------------------------------------------------
  ##   cT     ... parameter for transformation
  ##   x      ... argument
  ## -----------------------------------------------------------------------
  ## Return: antiderivative of inverse transformation at x.
  ## -----------------------------------------------------------------------

  ## Remove attributes (s.t. 'identical' works as expected).
  cT <- as.numeric(cT)

  if (identical(cT, 0))
    ## Case: T(x) = log(x)
    return (exp(x))

  if (identical(cT, -0.5))
    ## Case: T(x) = -1/sqrt(x)
    return (-1/x)
  
  if (identical(cT, -1))
    ## Case: T(x) = -1/x
    return (-log(-x))

  else
    ## Case: T(x) = sign(c) * x^c
    return (sign(cT) * cT/(cT+1) * (sign(cT) * x)^((cT+1)/cT))
}

## --------------------------------------------------------------------------

FTinv <- function(cT, x) {
  ## ------------------------------------------------------------------------
  ## Compute inverse of antiderivative of inverse transformation.
  ## ------------------------------------------------------------------------
  ##   cT     ... parameter for transformation
  ##   x      ... argument
  ## ------------------------------------------------------------------------
  ## Return: inverse of antiderivative of inverse transformation at x.
  ## ------------------------------------------------------------------------

  ## Remove attributes (s.t. 'identical' works as expected).
  cT <- as.numeric(cT)

  if (identical(cT, 0))
    ## Case: T(x) = log(x)
    return (log(x))

  if (identical(cT, -0.5))
    ## Case: T(x) = -1/sqrt(x)
    return (-1/x)

  if (identical(cT, -1))
    ## Case: T(x) = -1/x
    return (-exp(-x))

  else
    ## Case: T(x) = sign(c) * x^c
    return (sign(cT) * (sign(cT) * (cT+1)/cT * x)^(cT/(cT+1)))
}

## --------------------------------------------------------------------------
