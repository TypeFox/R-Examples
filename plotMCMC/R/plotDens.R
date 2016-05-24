plotDens <- function(mcmc, probs=c(0.025,0.975), points=FALSE, axes=TRUE, same.limits=FALSE,
                     between=list(x=axes,y=axes), div=1, log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL,
                     cex.main=1.2, cex.lab=1, cex.axis=0.8, cex.strip=0.8, col.strip="gray95", las=0, tck=0.5,
                     tick.number=5, lty.density=1, lwd.density=3, col.density="black", lty.median=2, lwd.median=1,
                     col.median="darkgray", lty.outer=3, lwd.outer=1, col.outer="darkgray", pch="|", cex.points=1,
                     col.points="black", plot=TRUE, ...)
{
  ## 1  Define functions
  panel.dens <- function(x, ...)
  {
    if(any(is.finite(x)) && var(x)>0)
      panel.densityplot(x, lty=lty.density, lwd=lwd.density, col.line=col.density, plot.points=points, pch=pch,
                        cex=cex.points, col=col.points, ...)
    else
      panel.densityplot(x, col="transparent", ...)
    panel.abline(v=quantile(x,probs=probs), lty=lty.outer, lwd=lwd.outer, col=col.outer)
    panel.abline(v=median(x), lty=lty.median, lwd=lwd.median, col=col.median)
  }

  ## 2  Parse args
  relation <- if(same.limits) "same" else "free"

  ## 3  Prepare data (rearrange, transform)
  if(is.null(dim(mcmc)))  # vector or mcmc(vector)
  {
    mcmc.name <- tail(as.character(substitute(mcmc)),1)
    mcmc <- matrix(mcmc, dimnames=list(NULL,mcmc.name))
  }
  mcmc <- if(log) log(mcmc/div,base=base) else mcmc/div
  mcmc <- as.data.frame(mcmc)
  n <- nrow(mcmc)
  p <- ncol(mcmc)
  x <- data.frame(Factor=ordered(rep(names(mcmc),each=n),names(mcmc)), Draw=rep(1:n,p), Value=as.vector(as.matrix(mcmc)))

  ## 4  Prepare plot (set pars, create list args)
  ocol <- trellis.par.get("strip.background")$col
  trellis.par.set(strip.background=list(col=col.strip))
  on.exit(trellis.par.set(strip.background=list(col=ocol)))
  mymain <- list(label=main, cex=cex.main)
  myxlab <- list(label=xlab, cex=cex.lab)
  myylab <- list(label=ylab, cex=cex.lab)
  myrot <- switch(as.character(las), "0"=0, "1"=0, "2"=90, "3"=90)
  myscales <- list(y=list(draw=FALSE, relation="free"),
                   x=list(draw=axes,relation=relation,cex=cex.axis,tck=tck,tick.number=tick.number,rot=myrot))
  mystrip <- list(cex=cex.strip)

  ## 5  Create trellis object
  graph <- densityplot(~Value|Factor, panel=panel.dens, data=x, as.table=TRUE, between=between,
                       main=mymain, xlab=myxlab, ylab=myylab, par.strip.text=mystrip, scales=myscales,
                       ...)
  if(!log)  # leave ylim alone if log-transformed
  {
    if(is.list(graph$y.limits))                                                 # set lower ylim to 0
      graph$y.limits <- lapply(graph$y.limits, function(y){y[1]<-0;return(y)})  # multi-panel plot
    else
      graph$y.limits[1] <- 0                                                    # single-panel plot
  }

  ## 6 Finish
  if(plot)
  {
    print(graph)
    invisible(x)
  }
  else
  {
    invisible(graph)
  }
}
