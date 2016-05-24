`plot.depratio` <-
  function (x, xlab = "Adjacent times", ylab = "Dependence ratio", 
            ylim=NULL, plot.ci = FALSE, legend = TRUE, ...) 
{
  if (class(x) != "depratio") 
    stop("plot only for depratio-objects")
  ntau <- dim(x$tau)
  if (plot.ci) {
    if (is.null(x$boot)) 
      stop("No bootstrap CI's: use boot.ci = TRUE")
    nam <- dimnames(x$tau)
    par(mfrow = c(ceiling(ntau[1]/2), ceiling(ntau[1]/2)))
    sapply(1:ntau[1], function(i) {
      plot(1:ntau[2], x$tau[i, ], xaxt = "n", 
           ylab = ylab, xlab = xlab,
           ylim = if(is.null(ylim))
           range(c(x$boot[[1]][i,],x$boot[[2]][i, ])) else(ylim), ...)
      segments(1:ntau[2], c(x$boot[[1]][i,]),1:ntau[2],c(x$boot[[2]][i, ]))
      if (ntau[1] > 1) 
        title(main = nam[[1]][i])
      abline(h = 1)
      axis(1, 1:ntau[2], labels = nam[[2]])
    })
    invisible(par(mfrow = c(1, 1)))
  }
  else {
    adds <- list(...)
    pch <- if (is.null(adds$pch)) 
      paste(seq(ntau[1]))
    else adds$pch
    lwd <- if (is.null(adds$lwd)) 
      1
    else adds$lwd
    cex <- if (is.null(adds$cex)) 
      1
    else adds$cex
    lty <- if (is.null(adds$lty)) 
      seq(min(ntau[1], 6))
    else adds$lty
    col <- if (is.null(adds$col)) 
      seq(min(ntau[1], 6))
    else adds$col
    ylim <- if (is.null(ylim)) 
      c(0.01, max(x$tau, na.rm = TRUE))
    else ylim
    matplot(seq(ntau[2]), t(x$tau), type = "b", 
            xaxt = "n", ylim = ylim, xlab = xlab, ylab = ylab, 
            ...)
      symbols(rep(seq(ntau[2]), ntau[1]), c(t(x$tau)), 
              circles = c(t(x$freq)), add = TRUE, inches = 0.2)
    axis(1, seq(ntau[2]), dimnames(x$tau)[[2]])
    abline(h = 1)
    if (legend)
      if (ntau[1] > 1) 
        legend(ncol(x$tau) - 1, ylim[2], dimnames(x$tau)[[1]], 
               pch = pch, lty = lty, col = col, lwd = lwd, merge = FALSE, 
               cex = cex)
  }
  invisible()
}
