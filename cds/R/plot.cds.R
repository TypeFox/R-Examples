#' Plot cds Objects
#' 
#' Plot method for \code{cds} objects
#' 
#' @param x An object of class \code{cds}.
#' @param which A numeric vector: a subset of \code{1:3} specifying the plots to produce.
#' @param type Passed to \code{\link{matplot}}.
#' @param lty Passed to \code{\link{matplot}}.
#' @param lwd Passed to \code{\link{matplot}}.
#' @param show.legend Logical; should a legend be added to the plot or not.
#' @param col Passed to \code{\link{matplot}}.
#' @param bty.legend Passed to \code{\link{legend}}.
#' @param intercept Logical indicating whether to plot the intercept.
#' @param scale Logical indicating whether an intercept should be included or not.
#' @param exp.factor Factor for expanding the x- and y-limits.
#' @param bubble.fact Passed to \code{\link{calc.wt.bubbles}} as argument \code{fact}.
#' @param cont.factor Continuity correction to apply in case one of the alpha's are
#' equal to zero.
#' @param add Logical; add to plot or not?
#' @param pch Plotting character to use.
#' @param \dots Additional arguments passed to \code{\link{points}}.
#' @keywords utility
#' @method plot cds
#' @export
plot.cds <- function(x, which = 1L:3L, type = "l", lty = 1, lwd = 2, show.legend = TRUE, 
                     col = colorspace::rainbow_hcl(nr), bty.legend = "n", 
                     intercept = ncol(x$alphamat) == 4, scale = FALSE, add = FALSE, 
                     exp.factor = 1.2, bubble.fact = 0.12, cont.factor = 0.01, pch = 15, ...){
  if (length(which) > 1) {
    oldpar <- par(ask = TRUE)
    on.exit(par(oldpar))
  }
  if (scale) intercept <- FALSE
  alphamat <- x$alphamat
  nr <- nrow(alphamat)
  if (!intercept & ncol(x$alphamat) == 4) alphamat <- x$alphamat[, -1] 
  if (1 %in% which) {
    xvals <- seq(from = 0, to = 1, length = 100)
    M <- ispline(xvals, tvec = c(0, 0.5, 1), intercept = intercept)
    coords <- tcrossprod(M, alphamat)
    if(scale) coords <- scale(coords, center = FALSE, scale = apply(coords, 2, max))
    matplot(xvals, coords, type = type, lty = lty, xlab = "Observed Scale", 
            ylab = "Optimal Scale", col = col, lwd = lwd, pch = pch, ...)
    if(show.legend) legend("bottomright", col = col, lwd = lwd, bty = bty.legend,
                           lty = lty, legend = rownames(alphamat))
  }
  
  if (2 %in% which){
    plot(x, which = 1)
    wts <- calc.wt.bubbles(dat = x$postrs, grp = x$grp, q = x$q, fact = bubble.fact)
    M <- ispline((1:x$q - 0.5)/x$q, tvec = c(0, 0.5, 1), intercept = intercept)
    bubcoords <- tcrossprod(M, alphamat) 
    symbols(rep((1:x$q - 0.5)/x$q, times = nr), as.numeric(bubcoords), 
            circles = as.numeric(t(wts)), add = TRUE, 
            fg = rep(col, each = x$q), bg = rep(col, each = x$q), inches = FALSE)
  }
  
  if (3 %in% which){
      if(any(alphamat == 0)) {
        alphamat <- alphamat + cont.factor
        warning("Continuity factor applied for alpha == 0")
      }
      curvpts <- cbind(alphamat[,"a2"]/alphamat[,"a1"], alphamat[,"a2"]/alphamat[,"a3"])
      curvpts <- log2(curvpts) - 1
      
  	if (add) {
  	points(curvpts, col = col, pch = pch, ...)
  	}
    if (!add) {
  	  curv.xlim <- range(c(curvpts[, 1], 0))
      curv.xlim <- exp.factor * c(-1, 1) * max(abs(curv.xlim))
      curv.ylim <- range(c(curvpts[, 2], 0))
      curv.ylim <- exp.factor * c(-1, 1) * max(abs(curv.ylim))
  	  plot(curvpts, xlim = curv.xlim, ylim = curv.ylim, xlab = expression(paste(log[2], 
            "(",alpha[2],"/",alpha[1],")", -1, sep="")), ylab = expression(paste(log[2],
            "(",alpha[2],"/",alpha[3],")", -1, sep="")), asp = 1, col = col, pch = pch, ...)
  	  abline(h = 0, col = "grey60", lwd = 2)
      abline(v = 0,  col = "grey60", lwd = 2)
  	}
  }
}
