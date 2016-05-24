quadset <- function(nquad=5, basisobj=NULL, breaks, returnMatrix=FALSE){

# last modified 8 May 2012 by Jim Ramsay

##
## 1.  Check nquad
##
  {
    if(nquad<5){
      warning("nquad must be at least 5;  increase to this minimum.")
      nquad <- 5
    }
    else {
      if((nquad%%2)!=1){
        warning("nquad must be an odd integer;  increased to enforce this.")
        nquad <- 1+2*floor(nquad/2)
      }
    }
  }
##
## 2.  check basisobj
##
  if(!is.null(basisobj) && !is.basis(basisobj))
    stop('basisobj is not a basis object.')
##
## 3.  check breaks
##
  if(missing(breaks) || length(breaks) == 0) {
    if(is.null(basisobj) || !is.basis(basisobj))
      stop("Either 'breaks' or 'basisobj' must be provided.")
#
    type <- basisobj$type
    if(type != 'bspline')
      stop(
        "'breaks' not supplied and 'basisobj' is not a spline basis.")
#
    rangeval <- basisobj$rangeval
    params   <- basisobj$params
    knots    <- c(rangeval[1], params, rangeval[2])
    breaks   <- unique(knots)
  }
##
## 4.  quadpts and quadwts
##
  nbreaks = length(breaks);
#
  db <- diff(breaks)
  nquad1 <- nquad-1
  nbreaks1 <- nbreaks-1
# 4.1.  First create quadpts as a matrix
  quadpts. <- array(NA, dim=c(nbreaks1, nquad) )
  quadpts.[, 1] <- breaks[-nbreaks]
  db. <- db/nquad1
  for(i in 2:nquad)
    quadpts.[, i] <- (quadpts.[, i-1]+db.)
# 4.2.  Now convert quadpts matrix to the desired vector
  quadpts <- as.vector(t(quadpts.))
# 4.3.  Similarly first create quadwts as a matrix
  quadwts. <- outer(c(1, 4, rep(c(2, 4), (nquad1-2)/2), 1),
                   db/(nquad1*3) )
# 4.4.  Now convert quadwts matrix to the desired vector
  quadwts <- as.vector(quadwts.)
  quadvals <- cbind(quadpts=quadpts, quadwts=quadwts)

  if(is.null(basisobj))return(quadvals)
#
  basisobj$quadvals <- quadvals
  values <- vector("list", 2)
  for( ivalue in 1:2){
    values[[ivalue]] <- eval.basis(quadpts, basisobj, ivalue-1, returnMatrix)
  }
  basisobj$values <- values
  basisobj
}

