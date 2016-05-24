plot.boxM <- 
  function(x, gplabel = NULL,
           pch=c(16, 15), 
           cex=c(2, 2.5), 
           col=c("blue", "red"),
           xlim,
           conf=0.95, method=2, bias.adj=TRUE, lwd=2, ...) {
    
    dets <- x$logDet
    ng <- length(dets)-1
    
    if (missing(xlim)) {
      xlim <- range(dets)
    }
    
    if (conf>0) {
      cov <- c(x$cov, list(pooled=x$pooled))
      n <- x$df + c(rep(1, ng), ng)
      CI <- logdetCI( cov, n=n, conf=conf, method=method, bias.adj=bias.adj )
      xlim[1] <- min(xlim[1], CI$lower)
      xlim[2] <- max(xlim[2], CI$upper)
    }
    
    dotchart(dets, xlab = "log determinant", xlim=xlim,  ...)
    if (conf>0) {
      arrows(CI$lower, 1:(ng+1), CI$upper, 1:(ng+1), lwd=lwd, angle=90, length=.075, code=3)
    }
    points(dets, 1:(ng+1),  
           cex=c(rep(cex[1], ng), cex[2]), 
           pch=c(rep(pch[1], ng), pch[2]),
           col=c(rep(col[1], ng), col[2]))
    
    if(!is.null(gplabel))
      text(par("usr")[1], ng + .5, gplabel, pos=4, cex=1.25)
    
  }
