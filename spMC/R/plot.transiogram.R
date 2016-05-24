plot.transiogram <-
function(x, ..., main, legend = FALSE, ci = NULL) {
  # Plot transition probabilities matrixes 1D
  #
  #          x transiogram object
  #        ... other args to pass to plot()
  #       main title string
  #     legend number of points per axes
  #         ci confidence intervall significance (e.g. `ci = 0.95`)

  n <- dim(x$Tmat)[1]
  if (!is.null(ci)) {
    if (ci >= 1 | ci <= 0) stop("\"ci\" must be a number between 0 and 1")
  }
  if (missing(main) || !is.character(main)) main <- "One-dimensional transiogram"
  
  ly <- matrix(c(rep(1, n), rep(0, n), 2:(n^2+1), rep(n^2+3, n), rep(n^2+4, n)), ncol = n, byrow = TRUE)
  yl0 <- c(rep(0, n + 4))
  yln <- c(rep(0, 2), rep(n^2 + 2, n), rep(0, 2))
  ly <- cbind(yln, yln, ly, yl0, yl0)
  widths <- c(3, 4/3.5, rep(30/n, n), 4/3.5, 3) / 4 
  heights <- c(0.75, 1/3.5, rep(7.5/n, n), 2.5/3.5, 1)
  if (!legend) {
    ly <- ly[-n^2-4, ]
    heights <- heights[-n^2-4]
  }
  ly <- layout(ly, widths = widths, heights = heights, respect = TRUE)

  par(mar=c(0.3, 0.1, 0.3, 0.1))
  plot.new()
  text(0.5, 0.5, main, cex = 2)

  nomi <- colnames(x$Tmat)

  UppBnd <- NULL #upper bound
  LowBnd <- NULL #lower bound
  if(!is.null(x$LOSE) & !is.null(ci)) {
    a <- .5 - ci * .5
    LogOdds <- .C('transLogOdds', mdim = as.integer(dim(x$LOSE)), 
                  empTR = as.double(x$Tmat), empTLO = double(prod(dim(x$LOSE))),
                  NAOK = TRUE, PACKAGE = "spMC")$empTLO
    UppBnd <- LogOdds - qnorm(a) * x$LOSE
    LowBnd <- LogOdds + qnorm(a) * x$LOSE
    UppBnd <- array(.C('LogOddstrans', mdim = as.integer(dim(x$LOSE)),
                       empTLO = as.double(UppBnd),  empTR = as.double(UppBnd),
                       NAOK = TRUE, PACKAGE = "spMC")$empTR, dim = dim(x$LOSE))
    LowBnd <- array(.C('LogOddstrans', mdim = as.integer(dim(x$LOSE)),
                       empTLO = as.double(LowBnd), empTR = as.double(LowBnd),
                       NAOK = TRUE, PACKAGE = "spMC")$empTR, dim = dim(x$LOSE))
  }

  par(mar = 2 * (rep(1.25, 4) + c(1, 1, 0, 0)) / n)
  xpos <- mean(range(na.omit(x$lags)))
  for(j in 1:n) {
    for(k in 1:n) {
      myseq <- na.omit(cbind(x$lags, x$Tmat[j, k, ]))
      if(dim(myseq)[1] != 0) {
        plot(myseq[, 1], myseq[, 2], ylab = "", xlab = "", ..., ylim = 0:1, axes = FALSE)
        if(!is.null(LowBnd)) {
          myseq <- na.omit(cbind(x$lags, LowBnd[j, k, ]))
          lines(myseq[, 1], myseq[, 2], lty = 2, ...)
        }
        if(!is.null(UppBnd)) {
          myseq <- na.omit(cbind(x$lags, UppBnd[j, k, ]))
          lines(myseq[, 1], myseq[, 2], lty = 2, ...)
        }
      }
      else {
        plot.new()
      }
      box()
      if (k == 1) axis(2)
      if (k == n) axis(4, .5, labels = nomi[j], tick = FALSE, font = 3)
      if (j == n) axis(1)
      if (j == 1) axis(3, xpos, labels = nomi[k], tick = FALSE, font = 3)
    }
  }
  par(mar = rep(.005, 4))
  plot.new()
  text(0.3, 0.5, "Transition probabilities", srt = 90, cex = 1.3)
  plot.new()
  text(0.5, 0.2, "Lags", cex = 1.3)

  if (legend) {
    args <- list(...)
    orgs <- list()
    orgs$col <- par("col")
    orgs$lty <- NA
    orgs$lwd <- NA
    orgs$pch <- NA
    orgs$pt.cex <- NA
    if (!is.null(args$col)) orgs$col <- args$col
    if (is.null(args$type)) {
      orgs$pch <- 1
      orgs$pt.cex <- 1
    }
    else {
      if (args$type == "p") {
        if (!is.null(args$pch)) orgs$pch <- args$pch else orgs$pch <- 1
        if (!is.null(args$cex)) orgs$pt.cex <- args$cex else orgs$pt.cex <- 1
      }
      if (args$type %in% c("b", "c", "o")) {
        if (!is.null(args$pch)) orgs$pch <- args$pch else orgs$pch <- 1
        if (!is.null(args$cex)) orgs$pt.cex <- args$cex else orgs$pt.cex <- 1
        if (!is.null(args$lty)) orgs$lty <- args$lty else orgs$lty <- 1
        if (!is.null(args$lwd)) orgs$lwd <- args$lwd else orgs$lwd <- 1
      }
      if (args$type %in% c("l", "h", "s", "S")) {
        if (!is.null(args$lty)) orgs$lty <- args$lty else orgs$lty <- 1
        if (!is.null(args$lwd)) orgs$lwd <- args$lwd else orgs$lwd <- 1
      }
    }
    par(mar=c(1.2, 0.1, 2.2, 0.1))
    plot.new()
    legend("center", legend = paste(x$type, "transiogram"),
           lty = orgs$lty, lwd = orgs$lwd, pch = orgs$pch, 
           pt.cex = orgs$pt.cex, col = orgs$col, bty = "n",
           cex = 1.5, box.lwd = 0, box.lty = 0, horiz = TRUE)
  }

}
