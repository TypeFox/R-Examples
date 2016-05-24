plot.asp=function(x, select=NULL, drv=0, bands=TRUE, level=0.95,
         grid=50, pages=0, plot=TRUE, ylim=NULL, xlab=NULL, ylab=NULL,
         scb.lwd=1, scb.lty="dotted", shade=FALSE, shade.col=grey(0.85),
         residuals=FALSE, residuals.col="steelblue", bayes=FALSE, rug=TRUE,...){
  if (residuals & drv==0) scbobject=scbM(object=x,div=1000,level=level,select=select,drv=drv,calc.stdev=T, bayes=bayes)
  else scbobject=scbM(object=x,div=1000,level=level,select=select,drv=drv,calc.stdev=F, bayes=bayes)
  plot(scbobject,bands=bands,grid=grid, plot=plot,pages=pages,select=select,ylim=ylim,xlab=xlab,ylab=ylab,scb.lwd=scb.lwd,scb.lty=scb.lty,shade=shade, shade.col=shade.col,residuals=residuals,residuals.col=residuals.col, bayes=bayes,rug=rug,...)
}
