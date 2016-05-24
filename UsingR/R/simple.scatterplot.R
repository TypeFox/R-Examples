##' make a simple scatter plot with histogram
##'
##' @param x data
##' @param y data
##' @param ... passed on
##' @return NULL
##'
##' @export
"simple.scatterplot" <-
  function(x,y,...) {
    def.par <- par(no.readonly = TRUE)# save default, for resetting...
    n<-length(x)
    xhist <- hist(x,sqrt(n), plot=FALSE)
    yhist <- hist(y, sqrt(n), plot=FALSE)
    top <- max(c(xhist$counts, yhist$counts))
    xrange <- c(min(x),max(x))
    yrange <- c(min(y),max(y))
    nf <- layout(matrix(c(2,0,1,3),2,2,TRUE), c(3,1), c(1,3), TRUE)
    layout.show(nf)
    
    par(mar=c(3,3,1,1))
    plot(x, y, xlim=xrange, ylim=yrange, xlab="x", ylab="y",...)
    abline(lm(y~x))
    par(mar=c(0,3,1,1))
    barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,col=gray(.95))
    par(mar=c(3,0,1,1))
    barplot(yhist$counts, axes=FALSE, xlim=c(0, top),
            space=0,col=gray(.95), horiz=TRUE)
    
    par(def.par)#- reset to default
  }
