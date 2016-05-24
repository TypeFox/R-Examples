knots <- function(Fn, ...)UseMethod('knots')

knots.basisfd <- function(Fn, interior = TRUE, ...) {
##
## 1.  object$type = 'bspline'?
##
  type <- Fn$type
  oName <- substring(deparse(substitute(Fn)), 1, 33)
  if(is.null(type))
    stop('is.null((', oName, ')$type);  must be "bspline"')
  if(type != 'bspline')
    stop('(', oName, ')$type) = ', type[1], ';  must be "bspline"')
##
## 2.  knots
##
  int <- Fn$params
  if(interior) return(int)
#
  nord <- norder(Fn)
  rng <- Fn$rangeval
  allKnots <- c(rep(rng[1], nord), int, rep(rng[2], nord))
  return(allKnots)
}

knots.fd <- function(Fn, interior=TRUE, ...){
  knots(Fn$basis, interior=interior, ...)
}
knots.fdSmooth <- function(Fn, interior=TRUE, ...){
  knots(Fn$fd, interior=interior, ...)
}
