`plotCEF` <-
function (typ=NA, fun=NA, dis=NA, a2=NA, c=NA, p1=NA, p2=p1, x=0:200/200, add=TRUE, xlim=c(0,1), ylim=c(0,1), plt.pt=TRUE, plt.ptann=TRUE, xlab=NA, ylab=NA, ...) {
  b <- CEF(typ=typ, fun=fun, dis=dis, a2=a2, c=c, p1=p1, p2=p2)
  if (!add) {
    plot(x, b(x), type="l", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...) 
    }
  else {
    points(x, b(x), type="l", ...)
    }
  if (plt.pt) {
    if (!is.na(p1)) {
      points(p1, p2, pch=16)
      if (plt.ptann) {
        text(p1+0.1, p2-0.03, paste(p1, ", ", p2)) 
        }
      }
    }
  }

