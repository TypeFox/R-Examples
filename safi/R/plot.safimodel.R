plot.safimodel <- function(x, ylim = NULL, ...) {
  s.m <- x
  if (class(s.m) != "safimodel") 
    stop("object class is not safimodel")
  normalized.coefficients <- s.m$normalized.coefficients
  split.points <- s.m$split.points
  d.f <- s.m$d.f
  variable.names <- s.m$variable.names
  plot.coef <- function(x, time, xlab = NULL, ylab = NULL, ...) {
    if (is.null(xlab)) 
      xlab <- "input function domain"
    if (is.null(ylab)) 
      ylab <- "normalized regression index"
    palette <- c(rev(colorRampPalette(c("white", "blue"))(100)), colorRampPalette(c("white", "red"))(100))
    scaled <- (x)/(2 * max(abs(x))) * 199 + 101
    col <- palette[scaled]
    barplot(x, width = diff(time), space = 0, col = col[1:length(x)], ylab = ylab, xlab = xlab, 
            xaxt = "n", ...)
    axis(1, at = seq(time[1], time[length(time)], length.out = 5), labels = seq(time[1], time[length(time)], 
                                                                                length.out = 5))
    abline(h = 0)
    box()
  }
  if (is.null(ylim)) {
    h <- max(sapply(normalized.coefficients, max))
    l <- min(sapply(normalized.coefficients, min))
    ylim <- c(min(l - 1/8 * abs(l), 0), max(h + 1/8 * abs(h), 0))
  }
  # run plot.coef for each pair of elements from norm.coef and split.points
  par(mfrow = c(1, d.f))
  invisible(mapply(plot.coef, normalized.coefficients, split.points, variable.names, MoreArgs = list(ylim = ylim, 
                                                                                                     ...)))
} 