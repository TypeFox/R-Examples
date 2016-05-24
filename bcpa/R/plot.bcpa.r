#' Plotting method for BCPA output
#' 
#' Plotting method for the output of a BCPA analysis with vertical break points, superimposed estimates of the partitioned mean and variance estimates and color-coded autocorrelation estimates.
#' 
#' @param x a \code{windowsweep} object, i.e. the output of the \code{\link{WindowSweep}} function.
#' @param type whether to plot smooth or flat bcpa output
#' @param threshold for smooth BCPA, this is the minimum number of windows that must have identified a given changepoint to be illustrated.
#' @param clusterwidth for flat BCPA, this is the temporal range within which change points are considered to be within the same cluster. 
#' @param col.cp,col.mean,col.sd color of the vertical change points, mean estimate, and prediction interval (mu +- sigma), respectively. 
#' @param pt.cex expansion coefficient for point sizes.
#' @param legend logical - whether to draw a legend or not. 
#' @param rho.where where to place the legend for the time-scale / auto-correlation.  Can be one of "nowhere", "top", "bottom", "left", "right", "topleft", "topright", "bottomright", "bottomleft"
#' @param mu.where where (and whether) to place the legend box for the mean - same options as for \code{rho.where}
#' @param ... other arguments to pass to the \code{plot} base function.
#' @author Eliezer Gurarie
#' @method plot bcpa
#' @seealso Plots output of the \code{\link{WindowSweep}} function. 
#' @examples 
#' if(!exists("Simp.ws"))
#' {
#'  data(Simp)
#'  Simp.ws <- WindowSweep(GetVT(Simp), "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
#' }
#' 
#' plot(Simp.ws)
#' # this actually provides basically the exact original changepoints
#' plot(Simp.ws, threshold=7)
#' # here's a flat analysis
#' plot(Simp.ws, type="flat", clusterwidth=3, legend=FALSE)

plot.bcpa <- 
  function (x, type = c("smooth", "flat")[1], threshold = 3, clusterwidth = 1, 
            col.cp = rgb(0.5, 0, 0.5, 0.5), pt.cex = 0.5, legend = TRUE, 
            rho.where = "topleft", mu.where = "nowhere", col.sd="red", col.mean="black",...) 
  {
    windowsweep <- x
    x <- windowsweep$x
    t.POSIX <- windowsweep$t.POSIX
    t <- windowsweep$t
    ws <- windowsweep$ws
    
    plot(t.POSIX, x, type = "n", ...)
    lines(t.POSIX, x, col = "grey")
    if (type == "smooth") {
      
      if("pp.smooth" %in% names(windowsweep)) 
        pp <- windowsweep$pp.smooth else pp <- PartitionParameters(windowsweep, type = "smooth")
      
      GoodBreaks <- ws$Break[ws$Model > 0]
      GoodBreaks <- as.data.frame(table(GoodBreaks))
      GoodBreaks <- data.frame(GoodBreaks, t.POSIX = t.POSIX[match(GoodBreaks[,1], windowsweep$t)])
      
      GoodBreaks[, 1] <- as.numeric(as.character(GoodBreaks[,1]))
      GoodBreaks <- GoodBreaks[GoodBreaks$Freq >= threshold,]
      
      abline(v = GoodBreaks[, 3], lwd = GoodBreaks[, 2]/threshold*2, col = col.cp)
      rho.scaled <- pp$rho.hat/max(pp$rho.hat, na.rm = 1)
      rho.int <- round(rho.scaled * 999 + 1)
      palette(topo.colors(1000))
      points(t.POSIX, x, pch = 21, col="darkgrey", bg = rho.int, cex = pt.cex, lwd=0.5)
      lines(t.POSIX, pp$mu.hat, lwd = 1.5, col.mean)
      lines(t.POSIX, pp$mu.hat + pp$s.hat, col = col.sd, lwd = 1.5)
      lines(t.POSIX, pp$mu.hat - pp$s.hat, col = col.sd, lwd = 1.5)
      rho.hat <- pp$rho.hat
    }
    
    if (type == "flat") {
      
      cpsummary <- ChangePointSummary(windowsweep, clusterwidth = clusterwidth)
      phases <- cpsummary$phases
      breaks <- cpsummary$breaks
      whichphase <- findInterval(t, phases$t0)
      rho.hat <- phases$rho.hat[whichphase]
      rho.int <- round(rho.hat/max(rho.hat, na.rm = TRUE) * 
                         999 + 1)
      palette(topo.colors(1000))
      points(t.POSIX, x, pch = 21, col="darkgrey", bg = rho.int, cex = pt.cex, lwd=0.5)
      
      closematch <- rep(NA, length=nrow(phases))
      for(i in 1:nrow(phases))
        closematch[i] <- which(abs(t-phases$t0[i]) == min(abs(t-phases$t0[i])))[1]
      
      phases$t0.POSIX <- t.POSIX[closematch]
      phases$t1.POSIX <- t.POSIX[c(closematch[-1], length(t))]
      
      t.mid <- (windowsweep$t[-1]+windowsweep$t[-length(windowsweep$t)])/2
      
      segments(phases$t0.POSIX, phases$mu.hat, phases$t1.POSIX, 
               phases$mu.hat, lwd = 1.5, col=col.mean)
      segments(phases$t0.POSIX, phases$mu.hat - phases$s.hat, 
               phases$t1.POSIX, phases$mu.hat - phases$s.hat, col = col.sd, 
               lwd = 1.5)
      segments(phases$t0.POSIX, phases$mu.hat + phases$s.hat, 
               phases$t1.POSIX, phases$mu.hat + phases$s.hat, col = col.sd, 
               lwd = 1.5)
      abline(v = phases$t0.POSIX[-1], lwd = breaks$size/max(breaks$size) * 4, col = col.cp)
    }
    if (legend) {
      legend.cols <- topo.colors(1000)[seq(1, 1000, length = 5)]
      legend.rhos <- signif(seq(0, max(rho.hat, na.rm = TRUE), 
                                length = 5), 2)
      if (rho.where != "nowhere") 
        legend(rho.where, bg = "white", legend = c(expression(hat(rho)), 
                                                   legend.rhos), pch = 19, ncol = 3, col = c(0, 
                                                                                             legend.cols), xjust = 0.5, yjust = 0.5)
      if (mu.where != "nowhere") 
        legend(mu.where, bg = "white", legend = c(expression(hat(mu)), 
                                                  expression(hat(mu) %+-% hat(sigma))), lty = 1, 
               lwd = 2:1, col = c("black", "red"), xjust = 0.5, 
               yjust = 0.5)
    }
    palette("default")
  }
