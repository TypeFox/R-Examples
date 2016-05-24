plot.cnv.intensities <-
function (x, my.colors = c("black", "red", "blue"), ylab = "Peak Intensity", xlab = c("individuals","Phenotype"),
    case.control, cex.leg = 0.8, dens.bw = "nrd0", dens.adjust = 1, ...)
{
    old.mfrow <- par("mfrow")
    old.mar <- par("mar")
    on.exit(par(mfrow = old.mfrow, mar = old.mar))
    mm <- matrix(c(2:1), nrow = 1, ncol = 2, byrow = TRUE)
    layout(mm, widths = c(1, 0.4))
    xx <- attr(x, "meanRatio")
    k <- attr(x, "k")
    mm <- attr(x, "means")
    par(mar = c(5.1, 0, 4.1, 2.1))
    if (is.null(attr(x, "batches"))){
      den <- density(xx, bw = dens.bw, adjust = dens.adjust)
      plot(den$y, den$x, type = "l", axes = FALSE, xlab = "density", ylab = "", col = "red")
      polygon(den$y, den$x, col = "red1")
      points(rep(0, k), mm, pch = 16, col = my.colors)
    } else {
      batches <- attr(x, "batches")
      bb <- sort(unique(batches))
      xmax <- -Inf
      y.coord <- NULL
      den <- list()
      for (i in 1:length(bb)){
        den[[i]] <- density(xx[batches==bb[i]], bw = dens.bw, adjust = dens.adjust)    
        den.i <- den[[i]]
        if (max(den.i$y) > xmax)
          xmax <- max(den.i$y)
        y.coord <- c(y.coord,den.i$x)  
      }
      plot(1, type = "n", axes = FALSE, xlab = "density", ylab = "", xlim=c(0,xmax), ylim = range(y.coord))
      for (i in 1:length(bb)){
        lines(den[[i]]$y, den[[i]]$x, lty = i, lwd = 2)
        points(rep(0, k), mm[i,], pch = 16, col = my.colors)      
      }
      legend("bottomright", c(paste("Batch:", bb)), bty = "n", cex = cex.leg, lwd = 2, lty = 1:length(bb))     
    }
    ll <- par("usr")[3:4]
    par(mar = c(5.1, 4.1, 4.1, 0))
    num.copies <- attr(x, "num.copies")
    x.ord <- sapply(x, function(x.i) which(attr(x,"num.copies")==x.i))
    if (missing(case.control)) {
        plot(xx, ylim = ll, yaxs="i", xlab = xlab[1], type = "n", ylab = ylab, ...)
        points(xx, col = my.colors[x.ord])
    } else {
        tt <- sort(unique(case.control))
        tt <- tt[!is.na(tt)]
        if (length(tt) == 1) {
          stop("case.control must have 2 differents values at least")
        }
        if (length(tt) > 2) {
          plot(case.control, xx, col = my.colors[x.ord], ylim = ll, yaxs="i", xlab = xlab[2], ylab = ylab, ...)
        }
        if (length(tt) == 2) {
          plot(xx, ylim = ll, yaxs="i", xlab = xlab[1], type = "n", ylab = ylab[1], ...)
          o <- case.control == tt[1]
          n <- sum(o)
          points(1:n, xx[o], col = my.colors[x.ord[o]], pch = 16)
          o <- case.control == tt[2]
          points((n + 1):(n + sum(o)), xx[o], col = my.colors[x.ord[o]], pch = 4)
          legend("bottomright", as.character(tt), pch = c(16, 4),
              title = "Case-control status", bty = "n", horiz = TRUE,
              cex = cex.leg)
        }
    }
    if (is.null(attr(x, "batches")))
      abline(h = mm, lty = 2, lwd = 2)
    else
      for (i in 1:length(bb))
        abline(h = mm[i,], lty = i, lwd = 2)
    legend("bottomleft", as.character(num.copies), col = my.colors[1:max(table(x))],
        horiz = TRUE, pch = 16, bty = "n", title = "copy number estimation",
        cex = cex.leg)
    if (is.null(attr(x, "batches")))
      if (!is.null(attr(x, "mixture")$P))
          legend("topleft", paste("Goodness-of-fit (p value):",
              round(attr(x, "mixture")$P, 5)), col = 0, horiz = TRUE, pch = 16,
              bty = "n", cex = cex.leg)
    invisible()
}

