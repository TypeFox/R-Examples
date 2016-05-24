TaylorSpline <- function(object, ...) {
  UseMethod('TaylorSpline')
}

TaylorSpline.dierckx <- function(object, ...) {
#  if(require(DierckxSpline)){
#    fdo <- dierckx2fd(object)
    fdo <- as.fd(object)
    return(TaylorSpline(fdo, ...))
#  }
#  else
#    stop('Requires library(DierckxSpline);  not installed.')
}

TaylorSpline.fd <- function(object, ...) {
##
## 1.  object$type = 'bspline'?
##
#  if(!require(fda))stop('fda package required.')
#
  oName <- substring(deparse(substitute(object)), 1, 33)
  type <- object$basis$type
  if(is.null(type))
    stop('is.null((', oName, ')$basis$type);  must be "bspline"')
  if(type != 'bspline')
    stop('(', oName, ')$basis$type) = ', type[1],
         ';  must be "bspline"')
##
## 2.  knots
##
  allKnots <- knots(object, interior=FALSE, ...)
  uniqKnots <- unique(allKnots)
  nUniq <- length(uniqKnots)
  nU1 <- (nUniq-1)
  midPts <- ((uniqKnots[-1]+uniqKnots[-nUniq])/2)
  nOrd <- norder(object)
##
## 3.  coef(object)
##
  coefObj <- as.array(coef(object))
  cdim <- dim(coefObj)
  cNames <- dimnames(coefObj)
  ndim <- length(cdim)
##
## 4.  switch(ndim, ...)
##
  colNames <- paste('b', 0:(nOrd-1), sep='')
  Dnames <- paste('D', 0:(nOrd-1), sep='')
  switch(ndim,
         '1'={
           Coef <- matrix(NA, nU1, nOrd, dimnames=list(
                                           NULL, colNames) )
           Deriv <- matrix(NA, nU1, nOrd, dimnames=list(
                                            NULL, Dnames) )
           for(i in 1:nOrd){
             bi <- eval.fd(midPts, object, i-1)
             Deriv[, i] <- bi
             Coef[, i] <- bi/factorial(i-1)
           }
         },
         '2'={
           Coef <- array(NA, c(nU1, nOrd, cdim[2]), dimnames=
                         list(NULL, colNames, NULL) )
           Deriv <- array(NA, c(nU1, nOrd, cdim[2]), dimnames=
                         list(NULL, Dnames, NULL) )
           if(!is.null(cNames) && !is.null(cNames[[2]])){
             dimnames(Coef)[[3]] <- cNames[[2]]
             dimnames(Deriv)[[3]] <- cNames[[2]]
           }
           for(i in 1:nOrd){
             bi <- eval.fd(midPts, object, i-1)
             bi. <- as.matrix(bi)
             Deriv[,i, ] <- bi.
             Coef[,i,] <- bi./factorial(i-1)
           }
         },
         '3'={
           Coef <- array(NA, c(nU1, nOrd, cdim[2:3]), dimnames=
                         list(NULL, colNames, NULL, NULL) )
           Deriv <- array(NA, c(nU1, nOrd, cdim[2:3]), dimnames=
                         list(NULL, Dnames, NULL, NULL) )
           if(!is.null(cNames)){
             if(!is.null(cNames[[2]])){
               dimnames(Coef)[[3]] <- cNames[[2]]
               dimnames(Deriv)[[3]] <- cNames[[2]]
             }
             if(!is.null(cNames[[3]])){
               dimnames(Coef)[[4]] <- cNames[[3]]
               dimnames(Deriv)[[4]] <- cNames[[3]]
             }
           }
           for(i in 1:nOrd){
             bi <- eval.fd(midPts, object, i-1)
             Deriv[,i, , ] <- bi
             Coef[,i, ,] <- bi/factorial(i-1)
           }
         },
         'other'=stop('coef(', oName, ' is neither a vector, nor ',
           'a matrix nor a 3-d array.')
         )
##
## 4.  Done
##
  Taylor <- list(knots=allKnots, midpoints=midPts, coef=Coef,
       deriv=Deriv)
  class(Taylor) <- 'Taylor'
  Taylor
}

#TaylorSpline.fd <- function(object, ...) {
#  if(depends(DierckxSpline)){
#    fdo <- fd2dierckx(object)
#    return(TaylorSpline(fdo, ...))
#  }
#  else
#    stop('Requires library(DierckxSpline;  not installed.')
#}

#TaylorSpline.list <- function(object, ...){
#  comps <- sapply(object, inherits, what)
#}

TaylorSpline.fdPar <- function(object, ...){
  TaylorSpline(object$fd, ...)
}

TaylorSpline.fdSmooth <- function(object, ...){
  TaylorSpline(object$fd, ...)
}

