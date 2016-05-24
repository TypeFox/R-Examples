plot_distr <- function(data, add = FALSE,
                       distr_col = adjustcolor("lightskyblue1", alpha.f = 0.4),
                       distr_border = "lightskyblue1",
                       bar_col = adjustcolor("azure3", alpha.f = 0.65),
                       bar_border = adjustcolor("azure3", alpha.f = 0.65), 
                       xlab = "", ylab = "", bars = TRUE, ...) {
  if (!add) {
    plot(x = data[ ,1], y = data[, 2], cex = 0, ylab = ylab, xlab = xlab, ...)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
         col = adjustcolor("grey", alpha.f = 0.30))
    axis(1, tck = 1, col.ticks = "white", labels = FALSE)
    axis(2, tck = 1, col.ticks = "white", labels = FALSE)
    if (bars)
      axis(3, at = data[, 1], labels = ifelse(as.integer(data[, 1]) == data[, 1], 
                                              as.integer(data[, 1]), round(data[, 1], 2)), 
           mgp = c(3, 0.35, 0), cex.axis = 0.85, tcl = -0.3, lwd.ticks = 1.2)
    
  }
  x <- c(data[1, 1], data[, 1], data[nrow(data), 1])
  y <- c(0, data[, 2], 0)
  polygon(x = x, y = y, col = distr_col, 
          border = NA)
  lines(x = x, y = y, col = distr_border, lwd = 2)
  if (bars) {
    bars <- calc_bars(data, w = (par("usr")[2] - par("usr")[1])/50)
    apply(bars, 1, function(x) 
      rect(x[1], x[2], x[3], x[4], col = bar_col, border = bar_border))
  }
}

#calculates coordinates of bars in our barcharts
calc_bars <- function(x, w = 1, names = "counts") {
  if (class(x)[1] == "matrix") {
    ytops <- x[ ,2]
    xs <- x[ ,1] 
  } else {
    if (names == "counts") {
      ytops <- as.vector(x)
      xs <- as.numeric(names(x))
    } else {
      ytops <- x
      xs <- 1L:length(x) - 1
    }
  }
  matrix(c(xs - 0.5*w, rep(0, length(ytops)), xs + 0.5*w, ytops), ncol = 4, byrow = F)
}
