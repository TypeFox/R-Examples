"boa.plot.gandr" <-
function(pname, bins = boa.par("gandr.bins"),
         alpha = boa.par("alpha"), win = boa.par("gandr.win"),
         annotate = boa.par("legend"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   drawn <- FALSE
   work <- boa.chain("work")
   work.support <- boa.chain("work.support")
   riter <- NULL
   for(i in names(work)) {
      if(is.element(pname, boa.pnames(work[[i]])))
         riter <- range(riter, boa.iter(work[[i]]))
   }
   x <- unique(round(seq(min(riter[1] + 49, riter[2]), riter[2],
                         length = bins)))
   R <- NULL
   Rq <- NULL
   for(i in x) {
      result <- boa.chain.gandr(work, work.support, alpha, pname, win, i)
      R <- c(R, result$csrf[1, 1])
      Rq <- c(Rq, result$csrf[1, 2])
   }
   idx <- is.finite(R)
   if(any(idx)) {
      drawn <- TRUE
      val <- boa.par("par")
      cex <- ifelse(is.null(val$cex), 1, val$cex)
      lwd <- ifelse(is.null(val$lwd), 1, val$lwd)
      x <- x[idx]
      R <- spline(x, R[idx])
      Rq <- spline(x, Rq[idx])
      ylim <- range(1, R$y, Rq$y)
      plot(R, xlab = "Last Iteration in Segment", ylab = "Shrink Factor",
           main = pname, ylim = ylim, type = "l", lwd = lwd)
      lines(Rq, lty = 2, lwd = lwd)
      abline(1, 0, lty = 3, lwd = lwd)
      usr <- par("usr")
      if(annotate)
         legend(x = usr[2], y = ylim[2], xjust = 1, yjust = 1,
                legend = c(paste(100 * (1 - alpha / 2), "%", sep=""), "Median"),
                lty = 2:1, bty = "n", cex = cex, lwd = lwd)
   }

   return(drawn)
}
