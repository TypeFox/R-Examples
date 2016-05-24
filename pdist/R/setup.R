setClass("pdist", representation(dist = "numeric",
                                 n = "numeric",
                                 p = "numeric"),
         S3methods = T)

#' extract parts of pdist
#'
#' @name [
#' @aliases [,pdist-method
#' @docType methods
#' @rdname extract-methods
setMethod("[", "pdist", function(x, i, j, ...) {
  if (j > x@p | j < 1) stop("index j out of bounds")
  if (i > x@n | i < 1) stop("index i out of bounds")
  if (missing(j)) j = 1:x@p
  x@dist[(i - 1) * x@p + j]
})

#' @method as.matrix pdist
#' @S3method as.matrix pdist
as.matrix.pdist = function(x, ...) matrix(x@dist, x@n, x@p, byrow=T)
