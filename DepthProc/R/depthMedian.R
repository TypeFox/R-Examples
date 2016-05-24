#' @title Depth median
#' @docType methods
#' @rdname depthMedian-methods
#' 
#' @param x object of class Depth or matrix.
#' @param ... arguments passed to \code{\link{depth}} function (e.g method).
#' 
#' @description
#' 
#' Return point with maximum depth function value. If multiple points have the same value, mean average of them will be returned.
#' 
#' @export
#' 
#' @examples
#' 
#' # depthMedian for matrix
#' x = matrix(rnorm(600), nc = 3)
#' depthMedian(x)
#' 
#' # depthMedian works with object of class Depth
#' dp = depth(x)
#' depthMedian(dp) 
#' 
setGeneric("depthMedian", function(x,...) standardGeneric("depthMedian"))

#' @rdname depthMedian-methods
#' @aliases depthMedian,matrix
#' @export
setMethod("depthMedian", "matrix", function(x,...)
{
  depths = depth(x,x,...)
  med = x[depths == max(depths),]
  if(ncol(x) != length(med)) med = colMeans(med)
  med
})

#' @rdname depthMedian-methods
#' @aliases depthMedian,data.frame
#' @export
setMethod("depthMedian", "data.frame", function(x,...)
{
  x = as.matrix(x)
  depthMedian(x, ...)
})

#' @rdname depthMedian-methods
#' @aliases depthMedian,Depth
#' @export
setMethod("depthMedian", "Depth", function(x)
{
  pos = which(x== max(x))
  med = x@u[pos,]
  if(ncol(x@u) != length(med)) med = colMeans(med)
  med  
})

