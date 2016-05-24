#' @encoding UTF-8
#' @title Detect Outliers
#' @description Perform an exploaratory test to detect \emph{outliers}. The quantity for \emph{min} reveals the minimum deviation from the mean, the integer in \emph{closest} highlights the position of the element. The quantity for \emph{max} is the maximum deviation from the mean, and the \code{farthest} integer is the position of such higher quantity.
#'
#' @param x A numeric object
#' @param index A numeric value to be considered in the computations
#'
#' @return Returns the minimum and maximum values, respectively preceded by their positions in the \code{vector}, \code{matrix} or \code{data.frame}.
#'
#' @references Dixon, W.J. (1950) Analysis of extreme values. \emph{Ann. Math. Stat.} \bold{21(4),} 488--506.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @seealso \link{winsorize} for diminishing the impact of outliers.
#'
#' @examples
#' outliers(x <- rnorm(20))
#'
#' #data frame:
#' age <- sample(1:100, 1000, rep=TRUE);
#' outliers(age)
#'
#' @export
`outliers` <-
  function(x, index=NULL) {
    if (is.data.frame(x)) {
      as.data.frame(sapply(x, outliers, index))
    } else if (is.matrix(x)) {
      apply(x, 2, outliers, index)
    } else if (is.list(x)) {
      lapply(x, outliers, index)
    } else if (is.vector(x)) {
      if (!is.null(index)) {
        if (!is.list(index)) {
          index <- list(index) # make sure index is a list
        }
        unsplit(outliers(split(x,index),index=NULL),index)
      } else {

        mu <- mean(x)
        dev <- abs(x - mu)
        closest <- which.min(dev)
        farthest <- which.max(dev)
        min <- dev[closest]
        max <- dev[ farthest]
        output <- data.frame(closest, min, farthest, max)
        return(output)
      }
    } else {
      cat("non-numeric argument to 'outlier'",class(x),"\n",sep="")
    }
  }
NULL
