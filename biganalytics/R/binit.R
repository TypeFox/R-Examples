#' Count elements appearing in bins of one or two variables
#' @description Provides preliminary counting functionality to 
#' eventually support graphical exploration or as an alternative to 
#' \code{table}.  Note the availability of \pkg{bigtabulate}.
#' @param x a \code{\link[bigmemory]{big.matrix}} or a \code{\link{matrix}}.
#' @param cols a vector of column indices or names of length 1 or 2.
#' @param breaks a number of bins to span the range from the maximum to the 
#' minimum, or a vector (1-variable case) or list of two vectors (2-variable 
#' case) where each vector is a triplet of min, max, and number of bins.
#' @return a list containing (a) a vector (1-variable case) or a matrix 
#' (2-variable case) of counts of the numbers of cases appearing in each of 
#' the bins, (b) description(s) of bin centers, and (c) description(s) of 
#' breaks between the bins.
#' @details The user may specify the number of bins to be used, of equal 
#' widths, spanning the range of the data (the default is 10 bins).  The user 
#' may also specify the range to be spanned along with the number of bins, in 
#' case a summary of a subrange of the data is desired.  Either univariate or 
#' bivariate counting is supported.
#'
#' The function uses left-closed intervals [a,b) except in the right-most bin,
#' where the interval is entirely closed.
#' @importFrom bigmemory attach.big.matrix is.big.matrix
#' @export
#' @examples
#' y <- matrix(rnorm(40), 20, 2)
#' y[1,1] <- NA
#' x <- as.big.matrix(y, type="double")
#' x[,]
#' binit(y, 1:2, list(c(-1,1,5), c(-1,1,2)))
#' binit(x, 1:2, list(c(-1,1,5), c(-1,1,2)))
#'
#' binit(y, 1:2)
#' binit(x, 1:2)
#'
#' binit(y, 1:2, 5)
#' binit(x, 1:2, 5)
#'
#' binit(y, 1)
#' binit(x, 1)
#'
#' x <- as.big.matrix(matrix(rnorm(400), 200, 2), type="double")
#' x[,1] <- x[,1] + 3
#' x.binit <- binit(x, 1:2)
#' filled.contour(round(x.binit$rowcenters,2), round(x.binit$colcenters,2),
#'                x.binit$counts, xlab="Variable 1",
#'                ylab="Variable 2")
binit <- function(x, cols, breaks=10)
  {
    if (!is.matrix(x) && !is.big.matrix(x))
      stop("Error in binit: x must be a matrix or a big.matrix.")
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    if (length(cols)<1 || length(cols)>2) {
      stop("Error in binit: only 1 or 2 columns is supported.")
    }
    if ( !is.list(breaks) && (length(breaks)==1 || length(breaks)==2) ) {
      if (is.numeric(breaks)) {
        usebreaks <- breaks
      } else { stop("Error in binit: breaks must be numeric.\n") }
      if (length(cols)==1) {
        if (!is.big.matrix(x)) {
          breaks <- c(min(x[,cols], na.rm=TRUE),
                      max(x[,cols], na.rm=TRUE), usebreaks[1])
        } else {
          breaks <- c(colmin(x, cols, na.rm=TRUE),
                      colmax(x, cols, na.rm=TRUE), usebreaks[1])
        }
      }
      if (length(cols)==2) {
        if (length(usebreaks)==1) usebreaks <- c(usebreaks, usebreaks)
        if (is.big.matrix(x)) {
          mins <- colmin(x, cols, na.rm=TRUE)
          maxs <- colmax(x, cols, na.rm=TRUE)
        } else {
          mins <- apply(x[,cols], 2, min, na.rm=TRUE)
          maxs <- apply(x[,cols], 2, max, na.rm=TRUE)
        }
        breaks <- list(c(mins[1], maxs[1], usebreaks[1]),
                       c(mins[2], maxs[2], usebreaks[2]))
      }
    }
    if (!is.list(breaks) && length(breaks)!=3)
      stop("Error in binit: incorrect specification of breaks.")
    if (is.list(breaks) & (length(breaks)!=2 ||
                           length(breaks[[1]])!=3 ||
                           length(breaks[[2]])!=3))
      stop("Error in binit: serious breaks problem.")

    if (is.list(breaks)) {
      if (is.big.matrix(x)) {
        ret = .Call("binit2BigMatrix", x@address,
          as.double(cols), as.double(breaks[[1]]), as.double(breaks[[2]]),
          PACKAGE="biganalytics")
      } else {
        if (is.integer(x)) {
          ret = .Call("binit2RIntMatrix", x,
            as.double(cols), as.double(breaks[[1]]), as.double(breaks[[2]]),
            PACKAGE="biganalytics")
        } else {
          ret = .Call("binit2RNumericMatrix", x,
            as.double(cols), as.double(breaks[[1]]), as.double(breaks[[2]]),
            PACKAGE="biganalytics")
        }
      }
      ret <- matrix(ret, breaks[[1]][3], breaks[[2]][3])
      rb <- seq(breaks[[1]][1], breaks[[1]][2], length.out=breaks[[1]][3]+1)
      rnames <- (rb[-length(rb)] + rb[-1]) / 2
      cb <- seq(breaks[[2]][1], breaks[[2]][2], length.out=breaks[[2]][3]+1)
      cnames <- (cb[-length(cb)] + cb[-1]) / 2
      ret <- list(counts=ret, rowcenters=rnames, colcenters=cnames,
                  rowbreaks=rb, colbreaks=cb)
    } else {
      if (is.big.matrix(x)) {
        ret = .Call("binit1BigMatrix", x@address,
          as.double(cols), as.double(breaks),
          PACKAGE="biganalytics")
      } else {
        if (is.integer(x)) {
          ret = .Call("binit1RIntMatrix", x,
            as.double(cols), as.double(breaks),
            PACKAGE="biganalytics")
        } else {
          ret = .Call("binit1RNumericMatrix", x,
            as.double(cols), as.double(breaks),
            PACKAGE="biganalytics")
        }
      }
      b <- seq(breaks[1], breaks[2], length.out=breaks[3]+1)
      rnames <- (b[-length(b)] + b[-1]) / 2
      ret <- list(counts=ret, centers=rnames, breaks=b)
    }

    return(ret)
  }


