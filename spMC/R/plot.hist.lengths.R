plot.hist.lengths <-
function (x, ...) {
  args <- list(...)
  nk <- length(x) - 2
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
  if (!is.null(args$main)) { title <- args$main } else {
    title <- "Histograms of stratum "
    if (x$log) title <- paste(title, "log-", sep = "")
    title <- paste(title, "lengths - Direction (", 
                   round(x$direction[1], 2), sep = "")
    if (length(x$direction) > 1) {
      for (i in 2:length(x$direction))
        title <- paste(title, round(x$direction[i], 2), sep = ", ")
    }
    title <- paste(title, ")", sep = "")
  }
  text(0.5, 0.5, labels = title, cex = 2)
  warn <- options("warn")$warn
  options("warn" = -1)
  for (i in 1:nk) {
    par(mar = c(5, 4, 3, 1))
    plot(x[[i]], ..., main = nomi[i])
  }
  options("warn" = warn)
}
