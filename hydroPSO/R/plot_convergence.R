# File plot_convergence.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2008-2011 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                        'plot_convergence'                                    # 
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates: 13-Ene-2012 ; 14-Nov-2012 ; 20-Nov-2012                             #        
################################################################################
 
                      
plot_convergence <- function(x,
                       verbose=TRUE,
                       col=c("black", "darkolivegreen"),
                       lty=c(1,3), 
                       lwd=c(2,2), 
                       main="Global Optimum & Normalized Swarm Radius vs Iteration Number", 
                       xlab="Iteration Number", 
                       ylab=c("Global Optimum", expression(delta[norm]) ), 
                       pch=c(15, 18), 
                       cex=1, 
                       cex.main=1.4, cex.axis=1.2, cex.lab=1.2, 
                       legend.pos="topright",                       
                       ...,                       
                       #### PNG options ### 
                       do.png=FALSE,
                       png.width=1500,
                       png.height=900,
                       png.res=90,
                       png.fname="ConvergenceMeasures.png"                  
                       
                       ) {
                         
                         
  # Checking that 'x' exists
  if ( missing(x) ) stop( "Missing argument: 'x'" )

  if (ncol(x)!=6) stop( paste("Invalid argument: 'ncol(x) != 6 (", ncol(x), "!=6)", sep="") )

  ############################ Getting the values ##############################
  # Iteration number
  iter <- x[,1]  
  # Global best
  gbest <- x[,2]  
  # Normalized swarm radius
  nradius <- x[,5]
  
  ############################     Plotting       ############################## 
  msg <- "[ Plotting convergence measures"
  if (do.png) msg <- paste(msg, " into '", basename(png.fname), sep="")
  msg <- paste(msg, "' ... ]", sep="")
  if (verbose) message(msg)    
 
  if (do.png) png(filename=png.fname, width=png.width, height=png.height, res=png.res)
  
  # Saving default plotting parameters
  old.par <- par(no.readonly=TRUE)
  if (!do.png) on.exit(par(old.par))
  
  oma <- c(2,   1, 2,   1)
  mar <- c(5, 4.5, 4, 4.5)+.1
  par(oma=oma, mar=mar)

  # Plotting the Gbest
  plot(iter, gbest, col=col[1], lty=lty[1], lwd=lwd[1], pch=pch[1], ylab=ylab[1],  
       xlab=xlab, cex=cex, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, 
       main=main, type="o", ylim=range(pretty(gbest)), panel.first=grid(), ...)
  ynticks <- par("yaxp")[3] + 1
  
  # Plotting the Normalized swarm radius as secondary Y-axis. 
  # Based on code from: http://rgraphics.limnology.wisc.edu/line.php
  par(new=TRUE, oma=oma)
  ylim <- range(nradius, na.rm=TRUE)
  plot(iter, nradius, ,xaxt="n", yaxt="n", ylim=ylim, 
       xlab="", ylab="", col=col[2], lty=lty[2], lwd=lwd[2], pch=pch[2], 
       cex=cex, type="o", ...)
  
  yat <- seq(ylim[1], ylim[2], length=ynticks)
  axis(4, col=col[2], col.axis=col[2], at= yat, labels=format(yat, scientific=TRUE, digits=4))
  mtext(ylab[2], side=4, line=3, cex=1.4*cex.lab, col=col[2])
  
  # Showing a legend
  legend(legend.pos, legend=ylab, lty=lty, col=col, pch=pch, lwd=lwd, bty="n", cex=1.4*cex.lab)
  
  if (do.png) dev.off()
  
}  # 'plot_convergence' END
