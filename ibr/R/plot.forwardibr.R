plot.forwardibr <- function(x,global=FALSE,...) {
  hasylab <- function(...) !all(is.na(pmatch(names(list(...)),"ylab")))
  hasxlab <- function(...) !all(is.na(pmatch(names(list(...)),"xlab")))
  if (is.null(dim(x))) return("Only one variable is retained") else {
    if (!global) {
      maxx <- apply(x,1,function(x) max(x[is.finite(x)]) )
      minn <- apply(x,1,function(x) min(x[is.finite(x)]) )
      x <- (x-matrix(minn,ncol=ncol(x),nrow=nrow(x)))/matrix(maxx,ncol=ncol(x),nrow=nrow(x))
    }
    x[!is.finite(x)] <- max(x[is.finite(x)])+0.1*diff(range(x[is.finite(x)]))
    if (hasylab(...)&hasxlab(...))
      
      image(1:ncol(x),1:nrow(x),t(x),...) else {
        if (hasxlab(...)) {
          image(1:ncol(x),1:nrow(x),t(x),ylab="Number of variables",...) }
        else {
          image(1:ncol(x),1:nrow(x),t(x),xlab="Variables",ylab="Number of variables",...) }
      }
    coordmin <- apply(x,1,which.min)
     text(coordmin,1:ncol(x),paste("V",1:ncol(x),sep=""))
    abline(h=seq(1.5,nrow(x)-0.5,by=1),col="grey10")
    abline(v=seq(1.5,ncol(x)-0.5,by=1),col="grey10")
  invisible(NULL)
  }
}
