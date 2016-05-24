plot.vario <-
  function (x, xlab="Lag distance", ylab=NULL, ylim=NULL, 
            xtype=c("mean.bin.dist", "bins"), rug=FALSE, ci=FALSE,
            pch=21, col.sig="black", col.nonsig="black", bg.sig="black", 
            bg.nonsig="white", alpha=0.05, ...) {
    
    xtypes=c("mean.bin.dist", "bins")
    xtype=match.arg(tolower(xtype), xtypes)
    
    switch (x$metric, 
            semivar = { yrange=range(x$vario, na.rm=TRUE); 
                        yname="Semivariance";
                        h=NA},
            cov = { yrange=range(x$vario, na.rm=TRUE); 
                    ylab="Covariance";
                    h=0},
            pearson = { yrange=c(-1,1); 
                        yname="Pearson correlation";
                        h=0},
            spearman = { yrange=c(-1,1); 
                         yname="Spearman correlation";
                         h=0},
            kendall = { yrange=c(-1,1); 
                        yname="Kendall correlation";
                        h=0},
            moran = { yrange=c(-1.1, 1.1); 
                      yname="Moran's I";
                      h=0},
            geary = { yrange=c(0, 2.2); 
                      yname="Geary's C";
                      h=0}
    )
    
    if (xtype=="mean.bin.dist")
      xvals=x$mean.bin.dist
    else
      xvals=x$bins
    if (is.null(ylab))
      ylab=yname
    if (is.null(ylim))
      ylim=yrange
    
    if (ci) {
      plot(xvals, x$vario, xlab=xlab, ylab=ylab, ylim=ylim, col=col.sig, type="n", ...)
      ci.vals=apply(x$rands, 2, quantile, c(alpha/2, 1-alpha/2))
      polygon (c(xvals, rev(xvals)), 
               c(ci.vals[1,], rev(ci.vals[2,])), col="lightgray")
      points(xvals, x$vario, xlab=xlab, ylab=ylab, ylim=ylim, col=col.sig, ...)
    }
    else {
      plot(xvals, x$vario, xlab=xlab, ylab=ylab, ylim=ylim, col=col.sig, ...)
    }
    
    if (!is.null(x$pvals)) {
      points(xvals[x$pvals >= alpha], x$vario[x$pvals >= alpha], bg=bg.nonsig, pch=pch, col=col.nonsig, ...)
      points(xvals[x$pvals < alpha], x$vario[x$pvals < alpha], bg=bg.sig, pch=pch, col=col.sig, ...)
    } 
    if (!is.na(h)) {
      if (x$is.centered)
        hv=h
      else
        hv=x$regional.mean
      abline(h=hv, lty=2)
    }
    
    if (rug)
      rug(jitter(rep(xvals, x$npoints)))
  }
