#' Utilities for big.matrix objects of package bigmemory
#' 
#' Extend the \pkg{bigmemory} package with various analytics.  In addition
#' to the more obvious summary statistics (see \code{\link{colmean}}, etc...),
#' \pkg{biganalytics} offers \code{\link{biglm.big.matrix}},
#' \code{\link{bigglm.big.matrix}}, \code{\link{bigkmeans}},
#' \code{\link{binit}}, and \code{apply} for 
#' \code{big.matrix} objects.  Some of the functions may be used with native 
#' \R objects, as well, providing gains in speed and memory-efficiency.
#' 
#' @name biganalytics-package
#' @aliases biganalytics biganalytics-package
#' @docType package
#' @import methods
#' @author Michael J. Kane and John W. Emerson.
#' 
#' Maintainers Michael J. Kane <bigmemoryauthors@gmail.com>
#' @references \url{http://www.bigmemory.org}
#' @keywords package
#' @examples
#' library(bigmemory)
#' 
#' x <- big.matrix(5, 2, type="integer", init=0,
#'                 dimnames=list(NULL, c("alpha", "beta")))
#' x
#' x[,]
#' x[,1] <- 1:5
#' x[,]
#' mean(x)
#' colmean(x)
#' summary(x)
#' apply(x, 1, mean)
NULL
