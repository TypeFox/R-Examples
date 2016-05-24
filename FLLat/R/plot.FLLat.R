plot.FLLat <- function(x,type=c("features","weights"),f.mar=c(5,3,4,2),
                       f.xlab="Probe",w.mar=c(3,5,0,2),
                       samp.names=1:ncol(x$Theta),hc.meth="complete",...) {

    if (!inherits(x,"FLLat")) {
        stop("'x' must be of class 'FLLat'")
    }
        
  type <- match.arg(type)
  J <- ncol(x$Beta)
  L <- nrow(x$Beta)
  S <- ncol(x$Theta)

  B.css <- colSums(x$Beta^2)
  B.po <- order(B.css,decreasing=T)

  if (type=="features") {
    T.ra <- rowMeans(x$Theta)
    y.range <- range(t(x$Beta)*sign(T.ra))

    n.plot.row <- ceiling(sqrt(J))
    par(mfcol=c(n.plot.row,ceiling(J/n.plot.row)),mar=f.mar)
    for (j in 1:J) {
      b.col <- ifelse(((x$Beta[,B.po[j]]*sign(T.ra[B.po[j]]))<0),3,2)
      plot(1:L,x$Beta[,B.po[j]]*sign(T.ra[B.po[j]]),main=paste("Feature",j),
           ylab="",type="h",ylim=y.range,xlab=f.xlab,col=b.col,...)
    }
  } else {
    hc.meths <- c("ward","single","complete","average","mcquitty","median",
                  "centroid")
    hc.meth <- match.arg(hc.meth,hc.meths)
    hm.col <- colorpanel(999,"blue","gray25","yellow")
    col.lim <- max(abs(x$Theta))
    preT <- x$Theta[B.po,,drop=F]
    newT <- preT*sign(rowMeans(preT))
    w.ylab <- paste("Feature",1:J)

    layout(c(1,2),heights=c(0.25,1))
    par(mar=c(0,w.mar[2],2,w.mar[4]))
    if (J==1) {
      x.ord <- 1:S
      frame()
    } else {
      hc <- hclust(as.dist(1-cor(newT)),method=hc.meth)
      x.ord <- hc$order
      newT <- newT[,x.ord]
      plot(as.dendrogram(hc),axes=F,xaxs="i",leaflab="none")
    }
    par(mar=w.mar)
    image(1:S,1:J,t(newT)[,J:1,drop=F],zlim=c(-col.lim,col.lim),axes=F,
          xlab="",ylab="",col=hm.col,...)
    axis(1,1:S,labels=samp.names[x.ord],tick=0,las=3)
    axis(2,1:J,labels=w.ylab[J:1],las=1,tick=0)
  }

}
