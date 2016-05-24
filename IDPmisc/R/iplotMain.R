### iplotMain.R

iplotMain <- function(main, cex.main, cex)
  ## internal functions for ipairs and ilagplot
  ## Version  2009-02-10
{
  ## writing title
  par(mar=rep(0,4), las=1, cex=cex)
  frame()
  plot.window(xlim=c(0,1),ylim=c(0,1))
  text(0.5, 0, labels=main, cex=cex.main, font=2, xpd=NA,
       adj=c(0.5,0))
} # PlotMain

