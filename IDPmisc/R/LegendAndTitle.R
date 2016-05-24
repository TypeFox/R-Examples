### LegendAndTitle.R

LegendAndTitle <- function(main,cex.main,border,colramp,zmax)
  ## internal functions for ipairs and ilagplot
  ## Version  2008-02-03
{
  ## writing title
  par(mar=rep(0,4), las=1)
  plot.new()
  plot.new()
  if(!is.null(main)) {
    plot(c(0,1),c(0,1),ty="n",bty="n", axes=FALSE, ann=FALSE)
    text(0.5, 0.5, labels=main, cex=cex.main, font=2)
  } else
  plot.new()

  ##-------------------
  ## drawing legend
  ## calculate colcut only when common zmax is defined
  if(is.null(zmax)) {
    warning(c("Separate color coding for each image\n"),call. = FALSE)
    border <- FALSE
    mycol <- c(par("bg"),colramp(100))
  }else {
    zmax <- ceiling(zmax)
    mycol <- c(par("bg"),colramp(zmax))
  }

  lev <- 0:length(mycol)
  par(mar=c(1,3,1,3), las=1)
  plot.new()
  plot.window(xlim=c(0, 1), ylim=range(lev,na.rm=TRUE), xaxs="i", yaxs="i")
  if(border)
    rect(0, lev[-length(lev)], 1, lev[-1], col=mycol)
  else
    rect(0, lev[-length(lev)], 1, lev[-1], col=mycol, border=mycol)
  box()
  if(is.null(zmax)) {
    axis(side=4, at=c(0.5,length(mycol)-0.5), labels=c("0","max"))
  }else {
    ap <- pretty(lev)
    axis(side=4, at=ap+0.5, labels=paste(ap))
  }
} # LegendAndTitle

