plotQuant <- function(mcmc, style="boxes", probs=c(0.025,0.975), axes=TRUE, names=NULL, ylim=NULL, yaxs="i", div=1,
                      log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, cex.axis=0.8, las=1, tck=-0.015,
                      tick.number=8, lty.median=1*(style!="bars"), lwd.median=1+1*(style!="boxes"), col.median="black",
                      lty.outer=1+2*(style=="lines"), lwd.outer=1, col.outer="black", pch=16, cex=0.8, col="black",
                      boxfill="darkgray", boxwex=0.5, staplewex=0.5, sfrac=0.005, mai=c(0.8,1,1,0.6),
                      mgp=list(bottom=c(2,0.4,0),left=c(3,0.6,0),top=c(0,0.6,0),right=c(0,0.6,0)), ...)
{
  ## 1  Set options
  opar <- par("mai", "mar", no.readonly=TRUE); on.exit(par(opar)); par(mai=mai)  # mai changes mar

  ## 2  Parse args
  if(is.mcmc.list(mcmc))
    mcmc <- as.mcmc(mcmc)
  if(is.mcmc(mcmc))
    mcmc <- as.data.frame(mcmc)
  style <- match.arg(style, c("bars","boxes","lines"))
  names <- if(is.null(names)) names(mcmc) else names
  probs <- c(min(probs), 0.25, 0.50, 0.75, max(probs))
  ellipsis <- as.list(substitute(list(...)))[-1]
  if(is.null(dim(mcmc))) stop("Argument 'mcmc' must contain more than one chain, arranged in columns.")

  ## 3  Prepare data (transform, quantiles)
  mcmc <- if(log) log(mcmc/div,base=base) else mcmc/div
  x <- seq(along=names(mcmc))
  y <- apply(mcmc, 2, quantile, probs=probs, na.rm=TRUE)  # support columns containing only NA values

  ## 4  Prepare plot (ylim)
  if(is.null(ylim))
    ylim <- if(!log && all(y[is.finite(y)]>0))
      c(0,1.04*max(y[is.finite(y)])) else range(y[is.finite(y)])+c(-0.04,0.04)*diff(range(y[is.finite(y)]))

  ## 4  Draw plot
  if(style == "bars")
  {
    plot(NA, xlim=c(0.5,ncol(mcmc)+0.5), ylim=ylim, yaxs=yaxs, axes=FALSE, ann=FALSE, ...)
    plotCI(x, y[3,], ui=y[5,], li=y[1,], pch=pch, cex=cex, col=col, gap=0, sfrac=sfrac, lty=lty.outer, lwd=lwd.outer,
           barcol=col.outer, add=TRUE, ...)
    lines(x, y[3,], lty=lty.median, lwd=lwd.median, col=col.median, ...)
  }
  else if(style == "boxes")
  {
    medlwd <- if(!is.null(ellipsis$medlwd)) ellipsis$medlwd else lwd.median
    whisklty <- if(!is.null(ellipsis$whisklty)) ellipsis$whisklty else lty.outer
    pars <- list(medlwd=medlwd, boxfill=boxfill, boxwex=boxwex, whisklty=whisklty, staplewex=staplewex)
    plot(NA, xlim=c(0.5,ncol(mcmc)+0.5), ylim=ylim, yaxs=yaxs, axes=FALSE, ann=FALSE, ...)
    bxp(list(stats=y,names=names), add=TRUE, axes=FALSE, pars=pars, ...)
  }
  else  # if(style == "lines")
  {
    matplot(x, t(y[c(1,5),]), ylim=ylim, yaxs=yaxs, type="l", lty=lty.outer, lwd=lwd.outer, col=col.outer, axes=FALSE,
            ann=FALSE, ...)
    lines(x, y[3,], lty=lty.median, lwd=lwd.median, col=col.median, ...)
  }

  ## 5  Annotate plot
  axes <- if(is.logical(axes)&&axes) 1:4 else axes
  yticks <- pretty(ylim,n=tick.number)[pretty(ylim,n=tick.number)>=ylim[1] & pretty(ylim,n=tick.number)<=ylim[2]]
  if(1 %in% axes) axis(1, cex.axis=cex.axis, las=las, tck=tck, mgp=mgp$bottom, at=x, labels=names, ...)
  if(2 %in% axes) axis(2, cex.axis=cex.axis, las=las, tck=tck, mgp=mgp$left, at=yticks, ...)
  if(3 %in% axes) axis(3, cex.axis=cex.axis, las=las, tck=tck, mgp=mgp$top, at=x, labels=names, ...)
  if(4 %in% axes) axis(4, cex.axis=cex.axis, las=las, tck=tck, mgp=mgp$right, at=yticks, ...)
  if(any(axes)) box()
  if(is.null(ellipsis$ann) || ellipsis$ann)
  {
    title(main=main, ...)
    title(xlab=xlab, mgp=mgp$bottom, ...)
    title(ylab=ylab, mgp=mgp$left, ...)
  }

  ## 6  Finish
  coords <- list(x=x, y=y)
  invisible(coords)
}
