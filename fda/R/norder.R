norder <- function(x, ...)UseMethod("norder")

norder.fd <- function(x, ...)norder.bspline(x$basis)

norder.basisfd <- function(x, ...)norder.bspline(x)

norder.default <- function(x, ...){
#
  xB. <- sapply(x, function(x){
    inherits(x, 'basisfd') || inherits(x, 'fd')
  } )
  xB <- which(xB.)
#
  {
    if(length(xB)<1)
      stop("input is not a 'basisfd' object and does not have ",
           "a 'basisfd' component.")
    else
      return(norder(x[[xB[1]]]))
  }
}

norder.bspline <- function(x, ...){
  if(!('type' %in% names(x))){
    xName <- substring(deparse(substitute(x)), 1, 33)
    stop('object ', xName, " does NOT have a 'type' component, ",
         "and therefore can NOT be a functional data object")
  }
  if(x$type != 'bspline'){
    xName <- substring(deparse(substitute(x)), 1, 33)
    stop('object ', xName, " is of type = ", x$type,
         ";  'norder' is only defined for type = 'bsline'")
  }
  with(x, nbasis - length(params))
}
