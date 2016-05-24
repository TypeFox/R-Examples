plot.density.lengths <-
function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", zero.line = TRUE, ...) {
    nk <- length(x) - 2
    if (is.null(main)) {
      main <- "Empirical densities of "
      if (!x$log) main <- paste(main, "stratum ", sep = "")
      if (x$log) main <- paste(main, "log-", sep = "")
      main <- paste(main, "lengths - Direction (", 
                     round(x$direction[1], 2), sep = "")
      if (length(x$direction) > 1) {
        for (i in 2:length(x$direction))
          main <- paste(main, round(x$direction[i], 2), sep = ", ")
      }
      main <- paste(main, ")", sep = "")
    }
    if (!is.null(xlab)) {
      if(length(xlab) < nk) xlab <- rep(xlab[1], nk)
    }
    else {
      xlab <- c()
      for (i in 1:nk) 
        xlab[i] <- paste("N =", x[[i]]$n, "  Log-Bandwidth =", formatC(x[[i]]$bw))
    }
    nomi <- names(x)[1:nk]
    mr <- ceiling(sqrt(nk))
    mc <- floor(mr - nk / mr)
    ly <- matrix(c(rep(1, mr), rep(0, mr), 2:(mr^2+1)), ncol = mr, byrow = TRUE)
    yl <- c(rep(0, mr + 2))
    ly <- cbind(yl, ly, yl)[, -(2+mr+1-mc)]
    widths <- c(rep(2 - 0.5 * mr, 1), rep(30/mr, mr), rep(4 - 0.5 * mr, 1)) / 4
    heights <- c(0.75, 1/3.5, rep(7.5/mr, mr-mc))
    ly <- layout(ly, widths = widths, heights = heights, respect = TRUE)
    mar <- sum(c(0, 1, 0, -0.5, -0.2, rep(0, 4), 0.2, rep(c(-0.1, rep(0, 3)), 2), 0)[1:mr])
    mar <- rep(mar, 4)
    par(mar = mar)
    plot.new()
    text(0.5, 0.5, labels = main, cex = 2)
    for (i in 1:nk) {
      par(mar = c(5, 4, 3, 2))
      plot.default(x[[i]], main = nomi[i], xlab = xlab[i], ylab = ylab, type = type, ...)
      if (zero.line) abline(h = 0, lwd = 0.1, col = "gray")
    }
}

