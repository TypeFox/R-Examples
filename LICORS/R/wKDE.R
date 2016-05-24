#' @title Weighted kernel density estimator (wKDE)
#' @aliases mv_wKDE
#'
#' @description 
#' \code{wKDE} gives a (weighted) kernel density estimate (KDE) for univariate data.
#'
#' If weights are not provided, all samples count equally.  It
#' evaluates on new data point by interpolation 
#' (using \code{\link[stats]{approx}}).
#' 
#' @param x data vector
#' @param eval.points points where the density should be evaluated. 
#' Default: \code{eval.points = x}.
#' @param weights vector of weights. Same length as \code{x}. 
#' Default: \code{weights=NULL} - equal weight for each sample.
#' @param kernel type of kernel. Default: \code{kernel='Gaussian'}. 
#' See \code{\link[stats]{density}} and \code{\link[locfit]{locfit.raw}} 
#' for additional options.
#' @param bw bandwidth. Either a character string indicating the 
#' method to use or a real number. 
#' Default: \code{bw="nrd0"}. Again see \code{\link[stats]{density}} for 
#' other options.
#' @return
#' A vector of length \code{length(eval.points)} (or \code{nrow(eval.points)}) 
#' with the probabilities of each point given the nonparametric fit on \code{x}.
#' @keywords distribution smooth
#' @export
#' @examples
#' ### Univariate example ###
#' xx = sort(c(rnorm(100, mean = 1), runif(100)))
#' plot(xx, wKDE(xx), type = "l")
#' yy = sort(runif(50, -1, 4)-1)
#' lines(yy, wKDE(xx, yy), col = 2)
#' 

wKDE <- function(x, eval.points = x, weights = NULL, 
                 kernel = "gaussian", bw = "nrd0") {
  x <- na.omit(x)
  if (!is.character(bw)){
    optimal.bw <- bw
  } else {
    if (is.null(weights)) {
      optimal.bw <- bw.nrd(x)
    } else {
      sel <- (weights > quantile(weights, 0.9))
      if (sum(sel) < 5) {
        optimal.bw <- bw.nrd(x)
      } else {
        optimal.bw <- bw.nrd(x[sel])
      }
    }
  }
  dens.object <- density(x, bw = optimal.bw, weights = weights, kernel = kernel)
  invisible(approx(dens.object$x, dens.object$y, eval.points)$y)
} 

#' @rdname wKDE
#' @description 
#' 
#' \code{mv_KDE} uses the \code{\link[locfit]{locfit.raw}} function in the 
#' \pkg{locfit} package to estimate KDEs for multivariate data. Note: Use this 
#' only for small dimensions, very slow otherwise.
#' @keywords distribution smooth
#' @export
#' @examples
#' ### Multivariate example ###
#' XX = matrix(rnorm(100), ncol = 2)
#' YY = matrix(runif(40), ncol = 2)
#' dens.object = mv_wKDE(XX)
#' 
#' plot(dens.object)
#' points(mv_wKDE(XX, YY), col = 2, ylab = "")
#
mv_wKDE <- function(x, eval.points = x, weights = NULL, kernel = "gaussian") {   
  if (is.null(weights)){
    weights <- 1
  }
  if (kernel == "gaussian"){
    kernel <- "gauss"
  }
  invisible(predict(locfit.raw(x, weights = weights, kern = kernel), eval.points))       
}


