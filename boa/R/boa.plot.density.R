"boa.plot.density" <-
function(lnames, pname, bandwidth = boa.par("bandwidth"),
                             window = boa.par("kernel"),
                             annotate = boa.par("legend"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   ipname <- list()
   xydensity <- list()
   xlim <- NULL
   ylim <- NULL
   key.names <- NULL
   lnames <- intersect(names(work), lnames)
   k <- 0
   for(i in lnames) {
      ipname[[i]] <- intersect(boa.pnames(work[[i]]), pname)
      for(j in ipname[[i]]) {
         k <- k + 1
         width <- bandwidth(work[[i]][, j])
         xydensity[[k]] <- density(work[[i]][, j], n = 100, width = width,
                                   window = window)
         xlim <- range(xlim, xydensity[[k]]$x)
         ylim <- range(ylim, xydensity[[k]]$y)
      }
      key.names <- c(key.names, substring(i, first = 1, last = 16))
   }
   drawn <- k > 0
   if(drawn) {
      val <- boa.par("par")
      cex <- ifelse(is.null(val$cex), 1, val$cex)
      lwd <- ifelse(is.null(val$lwd), 1, val$lwd)
      plot(xlim, ylim, xlab = pname, ylab = "Density", xlim = xlim,
           ylim = ylim, type = "n")
      k <- 0
      for(i in lnames) {
         for(j in ipname[[i]]) {
            k <- k + 1
            lines(xydensity[[k]], lty = k, lwd = lwd)
            parm <- work[[i]][, j]
         }
      }
      if(annotate)
         legend(x = xlim[2], y = ylim[2], xjust = 1, yjust = 1,
                legend = key.names, lty = 1:k, bty = "n",
                cex = cex, lwd = lwd)
   }

   return(drawn)
}
