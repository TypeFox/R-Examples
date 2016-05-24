`plotBounds` <-
function (a1=0, a0=1, add=TRUE, xlab=NA, ylab=NA, ...) {
  if (!add) {
    plot(.5, .5, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, col="transparent")
    }
  lines(c(a1, a1), c(0,1), ...) 
  lines(c(a0, a0), c(0,1), ...) 
  }

