#' Nearest neighborhoods for kernel smoothing
#' 
#' Nearest neighborhoods for the values of a continuous predictor. The result
#' is used for the conditional Kaplan-Meier estimator and other conditional
#' product limit estimators.
#' 
#' 
#' @param x Numeric vector -- typically the observations of a continuous random
#' variate.
#' @param bandwidth Controls the distance between neighbors in a neighborhood.
#' It can be a decimal, i.e.\ the bandwidth, or the string `"smooth"', in which
#' case \code{N^{-1/4}} is used, \code{N} being the sample size, or \code{NULL}
#' in which case the \code{\link{dpik}} function of the package KernSmooth is
#' used to find the optimal bandwidth.
#' @param kernel Only the rectangular kernel ("box") is implemented.
#' @return An object of class 'neighborhood'.  The value is a list that
#' includes the unique values of `x' (\code{values}) for which a neighborhood,
#' consisting of the nearest neighbors, is defined by the first neighbor
#' (\code{first.nbh}) of the usually very long vector \code{neighbors} and the
#' size of the neighborhood (\code{size.nbh}).
#' 
#' Further values are the arguments \code{bandwidth}, \code{kernel}, the total
#' sample size \code{n} and the number of unique values \code{nu}.
#' @author Thomas Gerds
#' @seealso \code{\link{dpik}}, \code{\link{prodlim}}
#' @references Stute, W. "Asymptotic Normality of Nearest Neighbor Regression
#' Function Estimates", \emph{The Annals of Statistics}, 1984,12,917--926.
#' @keywords smooth
#' @examples
#'
#' d <- SimSurv(20)
#' neighborhood(d$X2)
#' @export
"neighborhood" <- function(x,bandwidth=NULL,kernel="box"){
    if (any(is.na(x))) stop("Missing values in x")
    N <- length(x)
    if (N<2) stop("Not enough observations for kernel smoothing.")
    orderx <- order(x)
    values <- sort(unique(x))
    NU <- length(values)
    workx <- factor(x,labels=1:NU)
    tabu <- tabulate(workx)
    cumtabu <- cumsum(tabu)
    cumtabx <- rep(cumtabu,tabu)
    tabx <- rep(tabu,tabu)
    if (!length(bandwidth)){ ## need a bandwidth (dpik is from KernSmooth)
        ## require(KernSmooth)
        bandwidth <- KernSmooth::dpik(cumtabx/N,kernel="box")
    }
    else
        if (bandwidth=="smooth") bandwidth <- N^{-1/4}
    radius <- floor(bandwidth*N)
  
    nbh <- .C("neighborhood",
              first=integer(NU),
              size=integer(NU),
              as.integer(cumtabu),
              as.integer(cumtabx),
              as.integer(tabx),
              as.integer(radius),
              as.integer(NU),
              as.integer(N),
              PACKAGE="prodlim")
    nall <- sum(nbh$size)
    nbors <- .C("neighbors",
                first=nbh$first,
                size=nbh$size,
                as.integer(orderx),
                neighbors=integer(nall),
                as.integer(NU),
                PACKAGE="prodlim")$neighbors
  
  out <- list(values=values,
              first.nbh=nbh$first,
              size.nbh=nbh$size,
              neighbors=nbors,
              bandwidth=bandwidth,
              kernel=kernel,
              nu=NU,
              n=N)
    
  class(out) <- "neighborhood"
  out
}
