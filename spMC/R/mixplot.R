mixplot <-
function(x, main, legend = TRUE, ...) {
  # Plot transition probabilities matrixes 1D
  #
  #          x list of transiogram objects
  #       main title string
  #     legend logical for printing the legend
  #        ... other args to pass to plot()
  
  if (missing(main) || !is.character(main)) main <- "One-dimensional transiograms"
  if (is.transiogram(x)) x <- list(x)
  
  if (prod(sapply(x, is.transiogram)) == 0) stop("some elements of \"x\" are not transiogram objects")
  n <- sapply(x, function(n) dim(n$Tmat)[1])
  if (var(n) != 0) stop("transiogram objects must have the same number of categories")
  n <- n[1]
  nt <- length(x)

  ly <- matrix(c(rep(1, n), rep(0, n), 2:(n^2+1), rep(n^2+3, n), rep(n^2+4, n)), ncol = n, byrow = TRUE)
  yl0 <- c(rep(0, n + 4))
  yln <- c(rep(0, 2), rep(n^2 + 2, n), rep(0, 2))
  ly <- cbind(yln, yln, ly, yl0, yl0)
  widths <- c(3, 4/3.5, rep(30/n, n), 4/3.5, 3) / 4 
  heights <- c(0.75, 1/3.5, rep(7.5/n, n), 2.5/3.5, legend)
  ly <- layout(ly, widths = widths, heights = heights, respect = TRUE)
  
  par(mar=c(0.3, 0.1, 0.3, 0.1))
  plot.new()
  text(0.5, 0.5, main, cex = 2)
  nomi <- colnames(x[[1]]$Tmat)
  legends <- sapply(x, function(l) paste(l$type, "transiogram"))
  
  args <- list(...)
  if (!is.null(args$pch)) args$pch <- rep(args$pch, length.out = nt)
  if (!is.null(args$cex)) args$pt.cex <- rep(args$cex, length.out = nt)
  if (!is.null(args$lty)) args$lty <- rep(args$lty, length.out = nt)
  if (!is.null(args$lwd)) args$lwd <- rep(args$lwd, length.out = nt)
  orgs <- list()
  orgs$col <- par("col")
  orgs$lty <- rep(NA, length.out = nt)
  orgs$lwd <- rep(NA, length.out = nt)
  orgs$pch <- rep(NA, length.out = nt)
  orgs$pt.cex <- rep(NA, length.out = nt)
  if (!is.null(args$col)) orgs$col <- args$col
  orgs$col <- rep(orgs$col, length.out = nt)
  if (is.null(args$type)) {
    orgs$pch <- rep(1, length.out = nt)
    orgs$pt.cex <- rep(1, length.out = nt)
    orgs$type <- rep("p", length.out = nt)
  }
  else {
    orgs$type <- rep(args$type, length.out = nt)
    for (i in 1:nt) {
      if (orgs$type[i] == "p") {
        if (!is.null(args$pch)) orgs$pch[i] <- args$pch[i] else orgs$pch[i] <- 1
        if (!is.null(args$cex)) orgs$pt.cex[i] <- args$cex[i] else orgs$pt.cex[i] <- 1
      }
      if (orgs$type[i] %in% c("b", "c", "o")) {
        if (!is.null(args$pch)) orgs$pch[i] <- args$pch[i] else orgs$pch[i] <- 1
        if (!is.null(args$cex)) orgs$pt.cex[i] <- args$cex[i] else orgs$pt.cex[i] <- 1
        if (!is.null(args$lty)) orgs$lty[i] <- args$lty[i] else orgs$lty[i] <- 1
        if (!is.null(args$lwd)) orgs$lwd[i] <- args$lwd[i] else orgs$lwd[i] <- 1
      }
      if (orgs$type[i] %in% c("l", "h", "s", "S")) {
        if (!is.null(args$lty)) orgs$lty[i] <- args$lty[i] else orgs$lty[i] <- 1
        if (!is.null(args$lwd)) orgs$lwd[i] <- args$lwd[i] else orgs$lwd[i] <- 1
      }
    }
  }

  par(mar = 2 * (rep(1.25, 4) + c(1, 1, 0, 0)) / n)
  xlim <- range(as.vector(sapply(x, function(t) range(na.omit(t$lags)))))
  xpos <- mean(xlim)
  for(j in 1:n) {
    for(k in 1:n) {
      plot(xpos, .5, type = "n", ylab = "", xlab = "", xlim = xlim, ylim = 0:1, axes = FALSE)
      for (i in 1:nt) {
          myseq <- na.omit(cbind(x[[i]]$lags, x[[i]]$Tmat[j, k, ]))
        if (orgs$type[i] == "p") {
          points(myseq[, 1], myseq[, 2], type = "p", col = orgs$col[i],
                 cex = orgs$cex[i], pch = orgs$pch[i])      
        }
        if (orgs$type[i] %in% c("b", "c", "o")) {
          points(myseq[, 1], myseq[, 2], type = orgs$type[i],
               cex = orgs$cex[i], pch = orgs$pch[i], col = orgs$col[i],
               lty = orgs$lty[i], lwd = orgs$lwd[i])
        }
        if (orgs$type[i] %in% c("l", "s", "S", "h")) {
          points(myseq[, 1], myseq[, 2], type = orgs$type[i], 
                lty = orgs$lty[i], lwd = orgs$lwd[i], col = orgs$col[i])
        }
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
    par(mar=c(1.2, 0.1, 2.2, 0.1))
    plot.new()
    legend("center", legend = legends,
           lty = orgs$lty, lwd = orgs$lwd, pch = orgs$pch, 
           pt.cex = orgs$pt.cex, col = orgs$col, bty = "n",
           cex = 1.5, box.lwd = 0, box.lty = 0, horiz = TRUE)
  }
}
