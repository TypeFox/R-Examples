as.fd <- function(x, ...){
  UseMethod('as.fd')
}

as.fd.fdSmooth <- function(x, ...){
  x$fd
}

as.fd.dierckx <- function(x, ...){ 
# Translate an object of class dierckx to class fd
##
## 1.  check class
##
  objName <- deparse(substitute(x))
  {
    if(length(objName)>1)
      objName <- character(0)
    else
      if(nchar(objName)>33)
        objName <- substring(objName, 1, 33)
  }
  if(!inherits(x, 'dierckx')) 
    stop("'x' (", objName, ") is not of class 'dierckx', is ",
         class(x))
#    
  if(x$periodic)
    stop("object ", objName, " uses periodic B-splines.  ",
         "and as.fd is programmed to translate only ",
         "B-splines with coincident boundary knots.")
##
## 2.  Create a basis 
##
  rngval <- with(x, c(from, to))
  nKnots <- x$n
# length(dierckx$knots) = nest = estimated number of knots
# number actually used = dierckx$n  
#  knots <- object$knots[1:n]
  Knots <- knots(x, interior=FALSE)
  k <- x$k
  nOrder <- k+1
  breaks <- Knots[nOrder:(nKnots-k)]
#
  xlab <- x$xlab
#  if(!require(fda))
#    stop("library(fda) required for function 'dierckx2fd'",
#         ";  please install.") 
# 
  B.basis <- create.bspline.basis(rangeval=rngval, norder=nOrder,
                                breaks=breaks, names=xlab)
##
## 3.  Create fd object
##
  coef. <- coef(x)
#
  ylab <- x$ylab
  fdNms <- list(args=xlab, reps="reps 1", funs=ylab)
  fd(coef=coef., basisobj=B.basis, fdnames=fdNms)
}

as.fd.function <- function(x, ...){
# Translate an object of class splinefun to class fd
##
## 1.  check class
##
  objName <- deparse(substitute(x))
  {
    if(length(objName)>1)
      objName <- character(0)
    else
      if(nchar(objName)>33)
        objName <- substring(objName, 1, 33)
  }
  if(!inherits(x, 'function')) 
    stop("'x' (", objName, ") is not of class function")
#
  xenv <- environment(x)
  xz <- get('z', xenv) 
  if(is.null(xz))
    stop("NULL environment of 'x' (", objName,
         ");  therefore, it can NOT have been created by 'splinefun.'")
#  
  if(is.null(xz$method))
    stop("'x' (", objName, ") has a NULL 'method', and therefore",
         " can NOT have been created by 'splinefun.'")
# z$method:  1=periodic, 2=natural, 3=fmm (std B-Splines, I believe)   
#  if(xz$method!=3){
  if(!(xz$method %in% 2:3)){
    msg <- paste("x (", objName, ") ", sep='')
    msg2 <- {
      if(xz$method=="1")
        paste(msg, " uses periodic B-splines, and as.fd ",
              "is programmed\n    to translate only B-splines ",
              "with coincident boundary knots.", sep='')
      else
        paste(msg, "does not use B-splines as required ",
                      "for function 'as.fd'.")
    }
    stop(msg2)
  }
##
## 2.  Create a basis 
##
  Knots <- xz$x
  y.x <- xz$y
#  basis <- create.bspline.basis(breaks=Knots)
  basis <- create.bspline.basis(range(Knots), breaks=Knots)
  fd. <- fdPar(basis, lambda=0)
  nKn <- length(Knots) 
  nobs <- (2*nKn-1)
  x. <- seq(Knots[1], Knots[nKn], length=nobs) 
  smooth.basis(x., x(x.), fd.)$fd
}

as.fd.smooth.spline <- function(x, ...){
# Translate an object of class smooth.spline to class fd
##
## 1.  check class
##
  objName <- deparse(substitute(x))
  {
    if(length(objName)>1)
      objName <- character(0)
    else
      if(nchar(objName)>33)
        objName <- substring(objName, 1, 33)
  }
  if(!inherits(x, 'smooth.spline')) 
    stop("'x' (", objName, ") is not of class smooth.spline")
##
## 2.  Create a basis 
##
  Kn0 <- x$fit$knot
  x0 <- min(x$x)
  x1 <- max(x$x) 
  Knots <- (x0+(x1-x0)*Kn0[4:(length(Kn0)-3)])
# Don't use 'unique' in case 'x' has coincident interior knots.
#  basis <- create.bspline.basis(breaks=Knots)
  basis <- create.bspline.basis(range(Knots), breaks=Knots)
#
  fd(x$fit$coef, basis)
}

