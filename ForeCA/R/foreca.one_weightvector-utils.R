#' @title Plot, summary, and print methods for class 'foreca.one_weightvector'
#' @name foreca.one_weightvector-utils
#' @description
#' S3 methods for the one weightvector optimization in ForeCA 
#' (class \code{"foreca.one_weightvector"}).
#' @param x,object an object of class \code{"foreca.one_weightvector"}.
#' @param ... additional arguments passed to 
#' \code{\link[graphics]{plot}}, or \code{\link[base]{summary}}.
#' @param alpha significance level for testing white noise in 
#' \code{\link[stats]{Box.test}}; default: \eqn{0.05}.
#' @param lag integer; how many lags to test in \code{\link[stats]{Box.test}}; 
#' default: \eqn{10}.
#' @examples
#' # see examples in 'foreca.one_weightvector'
#'
NULL


#' @rdname foreca.one_weightvector-utils
#' @method summary foreca.one_weightvector
#' @description
#' \code{summary.foreca.one_weightvector} computes summary statistics.
#' @keywords manip
#' @export
#' 
summary.foreca.one_weightvector <- function(object, lag = 10, alpha = 0.05, ...) {
  .aux_pvalues <- function(xx) {
    Box.test(xx, lag = lag, type = "Ljung-Box")$p.value
  }
  pvals <- apply(object$score, 2, .aux_pvalues)
  
  out <- list(p.value = round(pvals, 4), 
              Omega = object$Omega, 
              alpha = alpha,
              lag = lag,
              weightvector = object$weightvector)
  
  out$selected <- which(out$p.value < alpha)
  return(out)
} 

#' @rdname foreca.one_weightvector-utils
#' @method plot foreca.one_weightvector
#' @description
#' \code{plot.foreca.one_weightvector} shows the results of an (iterative) 
#' algorithm that obtained the i-th optimal a weightvector 
#' \eqn{\mathbf{w}_i^*}.  It 
#' shows trace plots of the objective function and the weightvector, and a time series
#' plot of the transformed signal \eqn{y_t^*} along with its spectral density estimate 
#' \eqn{\widehat{f}_y(\omega_j)}.
#' @keywords manip hplot
#' @param main an overall title for the plot: see \code{\link[graphics]{title}}.
#' @param cex.lab size of the axes labels.
#' @export
#' 

plot.foreca.one_weightvector <- function(x, main = "", cex.lab = 1.1, ...) {
  
  object <- x
  temp.txt <- substitute(list(paste(hat(Omega) == w, "%")), 
                         list(w = round(object$Omega, 2)))
  
  object$best.f <- object$best.f + 1e-3  # add epsilon
  if (object$spectrum.control$smoothing && 
        requireNamespace("mgcv", quietly = TRUE)) {
    spec.tmp <- c(object$best.f)
    mod.tmp <- mgcv::gam(spec.tmp ~ s(seq_along(spec.tmp)), 
                         family = Gamma(link = "log"))
    # dispersion = 1 so we get an exponential distribution
    object$best.f.smoothed <- predict(mod.tmp, type = "response", dispersion = 1)
    object$best.f.smoothed <- normalize_mvspectrum(object$best.f.smoothed)
  } else {
    object$best.f.smoothed <- NULL
  }
  
  total.iter <- length(object$estimate$h.trace)
  
  par(mar = c(0, 4.5, 2, 0.5))
  layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE), heights = c(2, 3))

  plot(seq_len(total.iter) - 1, object$estimate$h.trace, 
       type = "l", lwd = 2, ylab = "", xlab = "", axes = FALSE, main = main, ...)
  axis(2)
  box()
  grid()
  points(seq_len(total.iter) - 1, object$estimate$h.trace, pch = 19)
  
  mtext(expression(paste("h(w|", hat(f)[U](omega[j]), ")")), side = 2, line = 2, 
        cex = cex.lab)

  plot(object$score, ylab = "", xlab = "", axes = FALSE, type = "l")
  box()
  grid()
  axis(2)
  
  par(mar = c(3.5, 4.5, 1, 0.5))
  matplot(seq_len(total.iter) - 1, object$estimate$weightvector.trace, type = "l", 
          ylab = "", xlab = "", axes = FALSE, lwd = 2, ...)
  axis(2)
  axis(1, at = seq_len(total.iter) - 1)
  mtext("weights", side = 2, line = 2, cex = cex.lab)
  mtext("Iteration", side = 1, line = 2.5, cex = cex.lab)
  box()
  grid()
  abline(h = 0)
  matpoints(seq_len(total.iter) - 1, object$estimate$weightvector.trace, 
            pch = 19, cex = cex.lab)
  
  plot(seq(0, 0.5, length = length(object$best.f) + 1)[-1], object$best.f, type = "l", 
       ylab = "", xlab = "", log = "y")
  lines(seq(0, 0.5, length = length(object$best.f) + 1)[-1], object$best.f.smoothed, 
        col = 4, lwd = 2)
  abline(h = 0.5, col = 2, lty = 2, lwd = 2)
  box()
  grid()
  mtext(expression(paste("Frequency / 2", pi)), side = 1, line = 2.5, cex = cex.lab)
  mtext(expression(paste(hat(f)(omega[j]), " (log scale)")), side = 2, 
        line = 2, cex = cex.lab)
  mtext(temp.txt, side = 3, adj = 1, line = -2, cex = cex.lab * 1.1)
} 

