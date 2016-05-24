#' @title Plot, summary, and print methods for class 'foreca'
#' @name foreca-utils
#' @description
#' A collection of S3 methods for estimated ForeCA results 
#' (class \code{"foreca"}).
#' @param x,object an object of class \code{"foreca"}.
#' @param ... additional arguments passed to 
#' \code{\link[stats]{biplot.princomp}}, \code{\link[stats]{biplot.default}},
#' \code{\link[graphics]{plot}}, or \code{\link[base]{summary}}.
#' @param alpha significance level for testing white noise in 
#' \code{\link[stats]{Box.test}}; default: \eqn{0.05}.
#' @param lag integer; how many lags to test in \code{\link[stats]{Box.test}}; 
#' default: \eqn{10}.
#' @examples
#' # see examples in 'foreca'
#'
NULL


#' @rdname foreca-utils
#' @method summary foreca
#' @description
#' \code{summary.foreca} computes summary statistics.
#' @keywords manip
#' @export
#' 
summary.foreca <- function(object, lag = 10, alpha = 0.05, ...) {
  .aux_pvalues <- function(series) {
    Box.test(series, lag = lag, type = "Ljung-Box")$p.value
  }
  
  pvals <- apply(object$scores, 2, .aux_pvalues)
  pvals.orig <- apply(object$series, 2, .aux_pvalues)
  
  out <- list(p.value = round(pvals, 4), 
              p.value.orig = round(pvals.orig, 4), 
              Omega = object$Omega,
              Omega.univ = object$Omega.univ,
              Omega.orig = 
                Omega(object$series, spectrum.control = object$spectrum.control,
                      entropy.control = object$entropy.control),
              alpha = alpha,
              lag = lag)
  
  out$selected <- which(out$p.value < alpha)
  return(out)
} 


#' @rdname foreca-utils
#' @method print foreca
#' @description
#' \code{print.foreca} prints a human-readable summary in the console.
#' 
#' @keywords manip
#' @export
#' 
print.foreca <- function(x, ...) {
  
  object <- x
  SO <- summary(object)
  
  cat("ForeCA found the top ", ncol(object$scores), " ForeCs of '",
      object$series.name, "' (", length(SO$Omega.orig), " time series).\n", sep ="")
  cat("Out of the top ", ncol(object$scores), " ForeCs, ", 
      ncol(object$scores) - length(SO$selected), 
      " are white noise.\n\n", sep = "")
  
  max.Omega.orig <- max(SO$Omega.orig)
  max.Omega.foreca <- SO$Omega.univ[1]
  
  cat("Omega(ForeC 1) = ", round(max.Omega.foreca, 2), "%",
      " vs. maximum Omega(", object$series.name, ") = ", 
      round(max.Omega.orig, 2), "%.\n", sep = "")
  omega.incr.abs <- max.Omega.foreca - max.Omega.orig 
  omega.incr.rel <- (max.Omega.foreca / max.Omega.orig - 1) * 100
  
  if (omega.incr.abs < 0) {
    warning("ForeCA converged to local optimum solution.\n\t ",
            "Please run foreca() again and increase the number of ",
            "random starts to avoid local minima.\n",
            "See also the 'Warning' section under ?foreca.")
  }
  
  cat("This is an absolute increase of ", round(omega.incr.abs, 2), " percentage points ",
      "(relative: ", round(omega.incr.rel, 2), "%) in forecastability.\n\n", sep = "")
  
  cat(rep("*", 10), "\n")
  cat("Use plot(), biplot(), and summary() for more details.\n\n")
} 

#' @rdname foreca-utils
#' @method biplot foreca
#' @description
#' \code{biplot.foreca} shows a biplot of the ForeCA loadings
#' (wrapper around \code{\link[stats]{biplot.princomp}}).
#' @keywords hplot
#' @export
#'

biplot.foreca <- function(x, ...) {
  object.princomp <- x
  class(object.princomp) <- "princomp"
  stats::biplot(object.princomp, ...)
  abline(h = 0)
  abline(v = 0)
}


#' @rdname foreca-utils
#' @method plot foreca
#' @description
#' \code{plot.foreca} shows biplots, screeplots, and white noise tests.
#' @keywords manip hplot
#' @export
#' 

plot.foreca <- function(x, lag = 10, alpha = 0.05, ...) {
  object <- x
  n.comp <- ncol(object$scores)

  SO <- summary(object, lag = lag, alpha = alpha)
    
  par(mar = c(4, 4.5, 2.5, 2.5))
  layout(matrix(1:6, ncol = 3, byrow = TRUE))
  biplot(object)
  
  ylim.max <- max(SO$Omega.orig, object$Omega.univ) * 1.05
  
  # replace 'Series' with 'ForeC'
  names(object$Omega.univ) <- gsub("Series ", "ForeC", names(object$Omega.univ))
  barplot(as.vector(object$Omega.univ), main = "Forecastability", 
          names.arg = names(object$Omega.univ), 
          ylab = expression(hat(Omega)(f[t]) ~ " (in %)"), 
          ylim = c(0, ylim.max))
  
  abline(h = object$Omega.univ[1], lty = 2, col = 4)
  abline(h = object$Omega.univ[n.comp], lty = 3, col = 4)
  
  barplot(as.vector(SO$Omega.orig), main = "Forecastability", 
          names.arg = names(SO$Omega.orig), 
          ylab = expression(hat(Omega)(x[t]) ~ " (in %)"), 
          ylim = c(0, ylim.max))
  abline(h = object$Omega.univ[1], lty = 2, col = 4)
  abline(h = object$Omega.univ[n.comp], lty = 3, col = 4)
  
  if (ncol(object$scores) >= 4) {
    biplot(object, 3:4)
  } else {
    plot.new()
  }
  barplot(SO$p.value, ylim = c(0, max(alpha * 1.1, SO$p.value)), ylab = "", 
          main = "")
  mtext("p-value \n (H0: white noise)", side = 2, line = 2, cex = 0.8)
  title(paste(sum(SO$p.value > alpha), "white noise"))
  # legend('topleft', c(paste(sum(SO$p.value > alpha), ' white noise')))
  abline(h = alpha, lwd = 2, lty = 2, col = 4)
  
  barplot(SO$p.value.orig, ylim = c(0, max(alpha * 1.1, SO$p.value.orig)), 
          ylab = "", main = "")
  mtext("p-value \n (H0: white noise)", side = 2, line = 2, cex = 0.8)
  title(paste(sum(SO$p.value.orig > alpha), "white noise"))
  abline(h = alpha, lwd = 2, lty = 2, col = 4)
}
