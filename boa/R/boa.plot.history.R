"boa.plot.history" <-
function(lnames, pname, annotate = boa.par("legend"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   ipname <- list()
   ybar <- list()
   xlim <- NULL
   ylim <- NULL
   key.names <- NULL
   lnames <- intersect(names(work), lnames)
   k <- 0
   for(i in lnames) {
      ipname[[i]] <- intersect(boa.pnames(work[[i]]), pname)
      for(j in ipname[[i]]) {
         k <- k + 1
         parm <- work[[i]][, j]
         ybar[[k]] <- cumsum(parm) / seq(parm)
         ylim <- range(ylim, ybar[[k]])
      }
      xlim <- range(xlim, boa.iter(work[[i]]))
      key.names <- c(key.names, substring(i, first = 1, last = 16))
   }
   drawn <- k > 0
   if(drawn) {
      val <- boa.par("par")
      cex <- ifelse(is.null(val$cex), 1, val$cex)
      lwd <- ifelse(is.null(val$lwd), 1, val$lwd)
      plot(xlim, ylim, xlab = "Iteration", ylab = pname,
           xlim = xlim, ylim = ylim, type = "n")
      k <- 0
      for(i in lnames) {
         for(j in ipname[[i]]) {
            k <- k + 1
            lines(boa.iter(work[[i]]), ybar[[k]], lty = k, lwd = lwd)
         }
      }
      if(annotate)
         legend(x = xlim[2], y = ylim[2], xjust = 1, yjust = 1,
                legend = key.names, lty = 1:k, bty = "n",
                cex = cex, lwd = lwd)
   }

   return(drawn)
}
