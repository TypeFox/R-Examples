`plot.pca` <- function(x, pc.axes=NULL, pch=16, col=par("col"), cex=0.8, mar=c(4, 4, 1, 1), ...) {

  opar <- par(no.readonly=TRUE)
  par(pty="s", cex=cex, mar=mar)
  
  if(is.null(pc.axes)) {
    ##-- Overview plot has been requested
    par(mfrow=c(2, 2))
    pc.axes <- 1:3
    
    ## Variance per PC for axis label annotation
    p <- paste0("PC",1:3," (",  round((x$L[1:3]/sum(x$L)) * 100, 2),"%)")

    plot(x$z[,1],x$z[,2], type="p", pch=pch, xlab=p[1], ylab=p[2], col=col, ...)
    abline(h=0,col="gray",lty=2); abline(v=0,col="gray",lty=2)

    plot(x$z[,3],x$z[,2], type="p", pch=pch, xlab=p[3], ylab=p[2], col=col,...)
    abline(h=0,col="gray",lty=2); abline(v=0,col="gray",lty=2)

    plot(x$z[,1],x$z[,3], type="p", pch=pch, xlab=p[1], ylab=p[3], col=col,... )
    abline(h=0,col="gray",lty=2); abline(v=0,col="gray",lty=2)

    plot.pca.scree(x$L, ...)
    par(opar)  ##- Reset par for subsequent plot to this device
    
  } else {
    ##-- Score plot of an individual PC pair (e.g. 1 vs 2) has been requested
    if(length(pc.axes) != 2) { 
      stop("Input 'pc.axes' should be NULL (for overview plots) or numeric and of length two (for individual score plot)") 
    }

    ## Variance for axis labels
    p <- paste0("PC",pc.axes," (", round((x$L[pc.axes]/sum(x$L)) * 100, 2),"%)")

    par(cex=cex+0.2) ## Make axis legend text a little larger
    plot(x$z[,pc.axes[1]],x$z[,pc.axes[2]], type="p", pch=pch, xlab=p[1], ylab=p[2], col=col, ...)
    abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)

  }
  invisible(x$z[,pc.axes, drop=FALSE])
}

