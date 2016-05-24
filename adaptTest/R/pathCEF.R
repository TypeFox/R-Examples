`pathCEF` <-
function (typ=NA, fun=NA, dis=NA, p1=1:49/50, p2=p1, x=0:200/200, plt.pt=FALSE, plt.ptann=FALSE, xlab=NA, ylab=NA, ...) { 
  plotCEF(typ=typ, fun=fun, dis=dis, p1=p1[1], p2=p2[1], x=x, add=FALSE, plt.pt=plt.pt, plt.ptann=plt.ptann, xlab=xlab, ylab=ylab, ...)
  for (i in 2:length(p1)) {
    plotCEF(typ=typ, fun=fun, dis=dis, p1=p1[i], p2=p2[i], x=x, add=TRUE, plt.pt=plt.pt, plt.ptann=plt.ptann, ...)
    }
  }

