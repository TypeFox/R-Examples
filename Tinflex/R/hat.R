#############################################################################
##
##  Parameters for hat and squeeze in an interval.
##
#############################################################################

##  Remark:
##  Tangents and secants have the form
##
##    t(x) = a + b*(x-y)
##
##  where 
##    a ... intercept for hat in transformed scale
##    b ... slope for hat in transformed scale
##    y ... anchor point of linear function
##
##  If 'y' is either the left boundary or the right boundary.
##  If possible then 'y' is chosen as the boundary point where the transformed
##  obtains higher valuation.

## --------------------------------------------------------------------------
## Codes for types of interval (see paper):

Ia   <- -1
Ib   <-  1
IIb  <-  2
IIa  <- -2
IIIa <- -3
IIIb <-  3
IVa  <- -4    ## concave
IVb  <-  4    ## convex

## --------------------------------------------------------------------------

type.iv <- function(left, right) {
  ## ------------------------------------------------------------------------
  ## Estimate type of interval.
  ## ------------------------------------------------------------------------
  ##   left   ... data for boundary point to the left
  ##   right  ... data for boundary point to the right
  ## ------------------------------------------------------------------------
  ## Return: type of interval;
  ##         0 if type cannot be estimated.
  ## ------------------------------------------------------------------------

  ## Parameter 'c' for interval.
  cT <- as.numeric(left["c"])
  
  ## Case: Unbounded domain.
  if ( is.infinite(left["x"])) {
    if (is.TRUE(right["d2Tfx"] < 0 && right["dTfx"] > 0))
      return (IVa)
    else
      return (0)
  }

  if ( is.infinite(right["x"])) {
    if (is.TRUE(left["d2Tfx"] < 0 && left["dTfx"] < 0))
      return (IVa)
    else
      return (0)
  }

  ## Case: Interval where the density vanishes at boundary
  ##       (and thus the log-density is -Inf).
  if ((is.TRUE(cT>0)  && identical(as.numeric(left["Tfx"]),0))    ||
      (is.TRUE(cT<=0) && identical(as.numeric(left["Tfx"]),-Inf)) ) {
    if (is.TRUE(right["d2Tfx"] < 0 && right["dTfx"] > 0))
      return (IVa)
    if (is.TRUE(right["d2Tfx"] > 0 && right["dTfx"] > 0))
      return (IVb)
    else
      return (0)
  }

  if ((is.TRUE(cT>0)  && identical(as.numeric(right["Tfx"]),0))    ||
      (is.TRUE(cT<=0) && identical(as.numeric(right["Tfx"]),-Inf)) ) {
    if (is.TRUE(left["d2Tfx"] < 0 && left["dTfx"] < 0))
      return (IVa)
    if (is.TRUE(left["d2Tfx"] > 0 && left["dTfx"] < 0))
      return (IVb)
    else
      return (0)
  }

  ## Case: Domains where the density has a pole at boundary
  ##       (and thus the transformed density equals 0 for c<0).
  if (is.TRUE(cT<0)) {
    if ((identical(as.numeric(left["Tfx"]),0)  && right["d2Tfx"] > 0) ||
        (identical(as.numeric(right["Tfx"]),0) && left["d2Tfx"] > 0)  )
      return (IVb)
  }

  ## Compute slope of secant.
  R <- (right["Tfx"] - left["Tfx"]) / (right["x"] - left["x"])

  ## Check for all other possible cases.
  if (is.TRUE(left["dTfx"] >= R && right["dTfx"] >= R))
    return (Ia)

  if (is.TRUE(left["dTfx"] <= R && right["dTfx"] <= R))
    return (Ib)

  if (is.TRUE(left["d2Tfx"] < 0 && right["d2Tfx"] < 0))
    return (IVa)

  if (is.TRUE(left["d2Tfx"] > 0 && right["d2Tfx"] > 0))
    return (IVb)

  if (is.TRUE(left["d2Tfx"] < 0 || right["d2Tfx"] > 0)) {
    if (is.TRUE(left["dTfx"] >= R && R >= right["dTfx"]))
      return (IIa)

    if (is.TRUE(left["dTfx"] <= R && R <= right["dTfx"]))
      return (IIIa)
  }

  if (is.TRUE(left["d2Tfx"] > 0 || right["d2Tfx"] < 0)) {
    if (is.TRUE(left["dTfx"] >= R && R >= right["dTfx"]))
      return (IIb)

    if (is.TRUE(left["dTfx"] <= R && R <= right["dTfx"]))
      return (IIIb)
  }

  ## Cannot estimate type of interval.
  ## Do we need a warning? Probably not.
  ## > warning("cannot detect type of interval")

  return (0)
}

## --------------------------------------------------------------------------

hat.iv <- function(left, right, link) {
  ## ------------------------------------------------------------------------
  ## Compute hat and squeeze for a paricular interval.
  ## ------------------------------------------------------------------------
  ##   left   ... data for boundary point to the left and for entire interval
  ##   right  ... data for boundary point to the right
  ##   link   ... pointer to next interval to the right
  ## ------------------------------------------------------------------------
  ## Return: vector with parameter of same kind as 'left'.
  ## ------------------------------------------------------------------------

  ## Insert pointer to next interval to the right
  ## (next element in poor man's linked list).
  left["next"] <- link

  ## Check for interval of length 0.
  if (is.TRUE( identical(left["x"],right["x"]) )) {
    left[c("ht.a","ht.b","ht.y")] <- c(left["Tfx"], 0, left["x"])
    left[c("sq.a","sq.b","sq.y")] <- c(left["Tfx"], 0, left["x"])
    left["A.ht"] <- 0
    left["A.sq"] <- 0
    left["type"] <- 0
    return (left)
  }
  
  ## Get type of distribution.
  left["type"] <- type <- type.iv(left, right)

  ## Compute tangent at boundary points.
  tl <- c( left["Tfx"], left["dTfx"], left["x"])
  tr <- c( right["Tfx"], right["dTfx"], right["x"])
  
  ## Compute secant.
  R <- (right["Tfx"] - left["Tfx"]) / (right["x"] - left["x"])
  if (is.TRUE(left["Tfx"] >= right["Tfx"])) 
    sc <- c( left["Tfx"], R, left["x"])
  else
    sc <- c( right["Tfx"], R, right["x"])
  
  ## Case: unbounded domains.
  if (is.infinite(left["x"]) && identical(type, IVa)) {
    left[c("ht.a","ht.b","ht.y")] <- tr
    left[c("sq.a","sq.b","sq.y")] <- NA
  }
  else if (is.infinite(right["x"]) && identical(type, IVa)) {
    left[c("ht.a","ht.b","ht.y")] <- tl
    left[c("sq.a","sq.b","sq.y")] <- NA
  }

  ## Case: bounded domains.
  else if (identical(type, Ia)) {
    left[c("ht.a","ht.b","ht.y")] <- tl
    left[c("sq.a","sq.b","sq.y")] <- tr
  }
  else if (identical(type, Ib)) {
    left[c("ht.a","ht.b","ht.y")] <- tr
    left[c("sq.a","sq.b","sq.y")] <- tl
  }
  else if (identical(type, IIa)) {
    left[c("ht.a","ht.b","ht.y")] <- tl
    left[c("sq.a","sq.b","sq.y")] <- sc
  }
  else if (identical(type, IIb)) {
    left[c("ht.a","ht.b","ht.y")] <- tr
    left[c("sq.a","sq.b","sq.y")] <- sc
  }
  else if (identical(type, IIIa)) {
    left[c("ht.a","ht.b","ht.y")] <- sc
    left[c("sq.a","sq.b","sq.y")] <- tr
  }
  else if (identical(type, IIIb)) {
    left[c("ht.a","ht.b","ht.y")] <- sc
    left[c("sq.a","sq.b","sq.y")] <- tl
  }
  else if (identical(type, IVb)) {
    left[c("ht.a","ht.b","ht.y")] <- sc
    if (is.TRUE(left["Tfx"] > right["Tfx"]))
      left[c("sq.a","sq.b","sq.y")] <- tr
    else
      left[c("sq.a","sq.b","sq.y")] <- tl
  }
  else if (identical(type, IVa)) {
    if (is.TRUE(left["Tfx"] > right["Tfx"]))
      left[c("ht.a","ht.b","ht.y")] <- tl
    else
      left[c("ht.a","ht.b","ht.y")] <- tr
    left[c("sq.a","sq.b","sq.y")] <- sc
  }
    
  ## Compute area below hat.
  left["A.ht"] <- area(left["c"], left["ht.a"], left["ht.b"], left["ht.y"], left["x"], right["x"])

  ## Compute area below squeeze.
  A.sq <- area(left["c"], left["sq.a"], left["sq.b"], left["sq.y"], left["x"], right["x"])
  left["A.sq"] <- if (is.finite(A.sq)) {A.sq} else {0} 

  ## Return vector with parameters.
  return (left)
}

## --------------------------------------------------------------------------

area <- function(cT, a,b,y, from,to) {
  ## -----------------------------------------------------------------------
  ## Compute area underneath hat or squeeze in particular interval.
  ## ------------------------------------------------------------------------
  ##   cT       ... parameter for transformation
  ##   a, b, y  ... intercept, slope and anchor point of transformed hat or squeeze
  ##   from, to ... range of integration
  ## ------------------------------------------------------------------------
  ## Return: area;
  ##         0 in case of an error.
  ## ------------------------------------------------------------------------

  ## Compute area ...
  area <- do.area(cT,a,b,y,from,to)

  ## ... and check result.
  
  if (! is.finite(area))
    ## We return Inf in all cases where 'area' is not finite (e.g. NaN or NA)
    area <- Inf

  if (area < 0) {
    ## Area is strictly negative.
    ## Two possible reasons:
    ## 1. The correct value of 'area' is extremely small or 0.
    ##    The negative number might result from cancelation.
    ##    Then 'area' can be set to 0.
    ## 2. Due to sever round-off error, the result is numerical garbage.
    ##
    ## We observed (2) in our experiments. So we discard the result and
    ## return 'Inf' (which means that the interval must be split).

    area <- Inf
  }

  return (area)
}

## ..........................................................................

do.area <- function(cT, a,b,y, from,to) {
  ## ------------------------------------------------------------------------
  ## Perform computation of area underneath hat or squeeze.
  ## ------------------------------------------------------------------------
  ##   cT       ... parameter for transformation
  ##   a, b, y  ... intercept, slope and anchor point of transformed hat or squeeze
  ##   from, to ... range of integration
  ## ------------------------------------------------------------------------
  ## Return: area.
  ## ------------------------------------------------------------------------

  ## Remove attributes (s.t. 'identical' works as expected).
  cT <- as.numeric(cT)

  ## Test where the tangent is constructed:
  ##   s = +1 if tangent is constructed on lower boundary of interval;
  ##       -1 if tangent is constructed on upper boundary of interval.
  s <- if (is.TRUE((to-y)>(y-from))) 1 else -1

  ## Generally we have
  ##   area <- (FT(cT, a+b*(to-y)) - FT(cT, a+b*(from-y))) / b
  ## For numerical reasons we have to distinguish
  ## between different values of 'cT'.

  if (identical(cT, 0)) {
    ## Case: T(x)=log(x)

    z <- s * b*(to-from)
    if (is.TRUE(abs(z) > 1.e-6)) {
      area <- (exp(a+b*(to-y)) - exp(a+b*(from-y))) / b
    } else {
      ## We need approximation by Taylor polynomial to avoid
      ## severe round-off errors.
      area <- exp(a) * (to-from) * (1 + z/2 + z*z/6)
    }
    return (area)
  }
  
  ## else: c!=0

  ## The tangent to the transformed density must not vanish.
  ## Otherwise, we cannot construct the hat function.
  ## Thus we simply return 'Inf' for the area.
  if (!(is.TRUE(sign(cT)*(a+b*(from-y)) >= 0)) ||
      !(is.TRUE(sign(cT)*(a+b*(to-y))   >= 0)) ) {
    return (Inf)
  }

  ## Transform b.
  z <- s * b/a * (to-from)
  
  if (identical(cT, -0.5)) {
    ## Case: T(x) = -1/sqrt(x)

    if (is.TRUE(abs(z) > 0.5)) {
      area <- (-1/(a+b*(to-y)) + 1/(a+b*(from-y))) / b
    } else {
      area <- 1/(a*a) * (to-from) / (1 + z)
    }
    return (area)
  }

  if (identical(cT, -1)) {
    ## Case: T(x) = -1/x

    if (is.TRUE(abs(z) > 1.e-6)) {
      area <- (-log(-a-b*(to-y)) + log(-a-b*(from-y))) / b
    } else {
      ## Approximation by Taylor polynomial.
      area <- -1/a * (to-from) * (1 - z/2 + z*z/3)
    }
    return (area)
  }

  if (identical(cT, 1)) {
    ## Case: T(x) = x

    area <- 0.5 * a * (to-from) * (z+2)
    return (area)
  }

  ## Case: T(x) = sgn(c) * x^c
  ## For all other cases we only use a rough approximation in
  ## case of numerical errors.
  
  if (is.TRUE(abs(b)>1e-10)) {
    area <- (FT(cT, a+b*(to-y)) - FT(cT, a+b*(from-y))) / b
  } else {
    area <- Tinv(cT, a) * (to-from)
  }
  return (area)

}

## --------------------------------------------------------------------------
