"boa.plot.bandg" <-
function(bins = boa.par("gandr.bins"), win = boa.par("gandr.win"),
         annotate = boa.par("legend"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   drawn <- FALSE
   work <- boa.chain("work")
   work.support <- boa.chain("work.support")
   riter <- NULL
   for(i in names(work))  riter <- range(riter, boa.iter(work[[i]]))
   x <- unique(round(seq(min(riter[1] + 49, riter[2]), riter[2],
                         length = bins)))
   Rp <- NULL
   Rmax <- NULL
   for(i in x) {
      result <- boa.chain.gandr(work, work.support, 1, window = win, to = i)
      Rp <- c(Rp, result$mpsrf)
      Rmax <- c(Rmax, max(result$psrf))
   }
   idx <- is.finite(Rp)
   if(any(idx)) {
      drawn <- TRUE
      val <- boa.par("par")
      cex <- ifelse(is.null(val$cex), 1, val$cex)
      lwd <- ifelse(is.null(val$lwd), 1, val$lwd)
      x <- x[idx]
      Rp <- spline(x, Rp[idx])
      Rmax <- spline(x, Rmax[idx])
      ylim <- range(1, Rp$y, Rmax$y)
      plot(Rmax, xlab = "Last Iteration in Segment", ylab = "Shrink Factor",
           ylim = ylim, type = "l", lwd = lwd)
      lines(Rp, lty = 2, lwd = lwd)
      abline(1, 0, lty = 3, lwd = lwd)
      usr <- par("usr")
      if(annotate)
         legend(x = usr[2], y = ylim[2], xjust = 1, yjust = 1,
                legend = c("Rp", "Rmax"), lty = c(2, 1), bty = "n",
                cex = cex, lwd = lwd)
   }

   return(drawn)
}
