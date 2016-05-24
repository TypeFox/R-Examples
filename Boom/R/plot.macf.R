PlotMacf <- function (x, lag.max = 40, gap = 0.5, main = NULL, boxes = TRUE,
    xlab = "lag", ylab = "ACF", type = "h") {
  ## Plot the autocorrelation function for a multivariate object.
  ## Args:
  ##   x: A numeric matrix, three-way array, or list of numeric
  ##     vectors.  If a matrix or array then the first index
  ##     corresponds to observation number.
  ##   lag.max:  The number of lags to compute in the ACF.
  ##   gap:  The size of the gap between plots.
  ##   main:  The main title for the plot.
  ##   boxes:  Logical.  Should boxes be drawn around the plots.
  ##   xlab:  The label for the horizontal axis.
  ##   ylab:  The label for the vertical axis.
  ##   type:  The plot type to use for plotting the ACF's.
  ##
  is.odd <- function(x) (x - 1)%%2 == 0

  if (is.list(x)) {
    nx <- length(x)
    nr <- max(1, floor(sqrt(nx)))
    nc <- ceiling(nx/nr)
  }
  else if (is.array(x)) {
    d <- dim(x)
    if (is.null(d) | length(d) == 1) {
      acf(as.vector(x), lag.max = lag.max, main = main,
          xlab = xlab, ylab = ylab)
      return(invisible(NULL))
    }
    ndim <- length(d)
    if (ndim == 3) {
      nr <- d[2]
      nc <- d[3]
      nx <- nr * nc
    }
    else if (ndim == 2) {
      nx <- d[2]
      nr <- max(c(floor(nx/sqrt(nx)), 1))
      nc <- ceiling(nx/nr)
    }
  }
  oma <- rep(4, 4)
  if (!is.null(main))
    oma[3] <- 6
  if (!is.null(xlab))
    oma[1] <- 6
  if (!is.null(ylab))
    oma[2] <- 6
  opar <- par(mfrow = c(nr, nc), mar = rep(gap/2, 4), oma = oma)
  on.exit(par(opar))
  a <- array(0, dim = c(lag.max + 1, nr, nc))
  m <- 0
  for (j in 1:nr) {
    for (k in 1:nc) {
      m <- m + 1
      if (m > nx)
        break
      if (is.list(x)) {
        y <- x[[m]]
        indx <- (1:length(y))[!is.na(y)]
        y <- y[indx]
      }
      else if (ndim == 3)
        y <- x[, j, k]
      else if (ndim == 2)
        y <- x[, m]
      else {
        cat("error in plot.macf:  wrong dimension for x\n")
        return(invisible(NULL))
      }
      if (length(y) == 1)
        tmp2 <- 1
      else tmp2 <- acf(y, lag.max = min(lag.max, length(y)),
                       plot = FALSE)$acf
      if (length(tmp2) < lag.max + 1) {
        zeros <- rep(0, lag.max + 1 - length(tmp2))
        tmp2 <- c(tmp2, zeros)
      }
      a[, j, k] <- tmp2
      if (any(is.na(a[, j, k]))) {
        if (sd(y) == 0) {
          a[1, j, k] <- 1
          a[2:(lag.max + 1), j, k] <- 0
        }
      }
    }
  }
  m <- 0
  lg <- 0:lag.max
  for (j in 1:nr) {
    ylim <- range(a[, j, ])
    for (k in 1:nc) {
      m <- m + 1
      if (m > nx) {
        plot(lg, type = "n", axes = FALSE)
        axis(1)
      }
      else {
        plot(lg, a[, j, k], axes = FALSE, type = type, ylim = ylim)
        if (boxes) {
          box()
        }
        abline(h = 0)
        if (j == nr)
          axis(1, xpd = NA)
        if (k == 1 & is.odd(j))
          axis(2, xpd = NA)
        if (k == nc & !is.odd(j))
          axis(4, xpd = NA)
      }
    }
  }
  if (!is.null(main))
    mtext(main, 3, 3, TRUE, 0.5, cex = par("cex.main"), font = par("font.main"))
  if (!is.null(xlab))
    mtext(xlab, 1, 3, TRUE, 0.5)
  if (!is.null(ylab))
    mtext(ylab, 2, 3, TRUE, 0.5)
  return(invisible(NULL))
}
