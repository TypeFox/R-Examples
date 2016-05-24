PlotSection <- function(var, horiz, depth, toptitle = '', sizetit = 1, 
                        units = '', brks = NULL, cols = NULL, axelab = TRUE, 
                        intydep = 200, intxhoriz = 20, drawleg = TRUE) {
  #
  #  Input arguments 
  # ~~~~~~~~~~~~~~~~~
  #
  dims <- dim(var)
  if (length(dims) > 2) {
    stop("Only 2 dimensions expected for var : (lon,depth) or (lat,depth)")
  }
  if (dims[1] != length(horiz) | dims[2] != length(depth)) {
    if (dims[1] == length(depth) & dims[2] == length(horiz)) {
      var <- t(var)
      dims <- dim(var)
    } else {
      stop("Inconsistent var dimensions and longitudes/latitudes + depth")
    }
  }
  dhoriz <- horiz[2:dims[1]] - horiz[1:(dims[1] - 1)]
  wher <- which(dhoriz > (mean(dhoriz) + 5))
  if (length(wher) > 0) {
    horiz[(wher + 1):dims[1]] <- horiz[(wher + 1):dims[1]] - 360
  }
  horizb <- sort(horiz, index.return = TRUE)
  depthb <- sort(-abs(depth), index.return = TRUE)
  horizmin <- floor(min(horiz) / 10) * 10
  horizmax <- ceiling(max(horiz) / 10) * 10
  depmin <- min(depth)
  depmax <- max(depth)
  if (is.null(brks) == TRUE) {
    ll <- signif(min(var, na.rm = TRUE), 4)
    ul <- signif(max(var, na.rm = TRUE), 4)
    if (is.null(cols) == TRUE) {
      cols <- c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
                "white", "white", "yellow", "orange", "red", "saddlebrown")
    }
    nlev <- length(cols)
    brks <- signif(seq(ll, ul, (ul - ll) / nlev), 4)
  } else {
    if (is.null(cols) == TRUE) {
      nlev <- length(brks) - 1
      cols <- rainbow(nlev)
    } else {
      if (length(cols) != (length(brks) - 1)) {
        stop("Inconsistent colour levels / list of colours")
      }
    }
  }
  #
  #  Plotting the section
  # ~~~~~~~~~~~~~~~~~~
  #
  xmargin <- 0.5
  ymargin <- 0.5
  topmargin <- 3
  if (axelab) {
    ymargin <- ymargin + 2.5
    xmargin <- xmargin + 1.5
  }
  if (drawleg) {
    layout(matrix(1:2, ncol = 1, nrow = 2), heights = c(5, 1))
    xmargin <- max(xmargin - 1.8, 0)
  }
  if (toptitle == '') { 
    topmargin <- topmargin - 2.5 
  } 
  par(mar = c(xmargin, ymargin, topmargin, 0.5), cex = 1.4, 
      mgp = c(2.5, 0.5, 0), las = 1)
  image(horizb$x, depthb$x, array(0, dims), col = 'grey', breaks = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "", main = toptitle, 
        cex.main = 1.5 * sizetit)
  image(horizb$x, depthb$x, var[horizb$ix, depthb$ix], col = cols, 
        breaks = brks, axes = FALSE, xlab = "", ylab = "", add = TRUE)
  if (axelab) {
    minhoriz <- ceiling(round(min(horizb$x), 0) / 10) * 10
    maxhoriz <- floor(round(max(horizb$x), 0) / 10) * 10
    axis(1, at = seq(minhoriz, maxhoriz, intxhoriz), tck = -0.02)
    maxdepth <- floor(round(max(depthb$x), 0) / 10) * 10
    axis(2, at = seq(-8000, 0, intydep), tck = -0.015)
  }
  box()
  #
  #  Colorbar
  # ~~~~~~~~~~
  #
  if (drawleg) {
    par(mar = c(1.5, ymargin, 2.5, 0.5), mgp = c(1.5, 0.3, 0), las = 1, 
        cex = 1.2)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = units, ylab = '')
    box()
    axis(1, at = seq(0.5, length(brks) - 0.5, 1), labels = brks, cex.axis = 1)
  }
}
