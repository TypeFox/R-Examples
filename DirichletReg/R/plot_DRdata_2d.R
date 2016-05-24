plot_DRdata_2d <- function(y, rug, main, ylim, colr, lwd, lty){
  colr <- if(is.null(colr)) "black" else colr
  y.dens <- density(y, from=0, to=1)
  ylims <- c(ifelse(is.null(ylim[1]), 0, ylim[1]),
             ifelse(is.null(ylim[2]), max(y.dens$y), ylim[2]))
  plot(y.dens, type="l", ylim=ylims, main=main, col=colr, lwd=lwd, lty=lty)
  if(rug) rug(y)
}
