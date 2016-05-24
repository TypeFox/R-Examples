`compHclust.heatmap` <-
function(x,xhc,gi,d.title="Cluster Dendrogram",hm.lab=TRUE,hm.lab.cex=1,
         d.ht=0.25,gi.width=0.5,d.mar=c(0,4,4,2),hm.mar=c(5,4,2,2)) {
  if (class(xhc)!="hclust") {
    stop("'xhc' must be of class 'hclust'")
  }
  if (ncol(x)!=length(xhc$order)) {
    stop("'xhc' must be a clustering of the columns of 'x'")
  }
  if (length(gi)!=nrow(x)) {
    stop("Length of 'gi' must equal number of rows of 'x'")
  }
  layout(cbind(c(1,2),c(3,4)),heights=c(d.ht,1),widths=c(1,gi.width))
  par(mar=d.mar)
  plot(as.dendrogram(xhc),axes=FALSE,xaxs="i",leaflab="none",main=d.title)
  par(mar=hm.mar)
  image(1:ncol(x),1:nrow(x),t(x[,xhc$order])[,nrow(x):1],axes=FALSE,xlab="",
        ylab="")
  if (hm.lab) {
    axis(3,1:ncol(x),labels=xhc$order,tick=0,cex.axis=hm.lab.cex)
  }
  par(mar=d.mar)
  frame()
  par(mar=hm.mar,yaxs="i")
  barplot(rev(gi),horiz=TRUE,main="Gene Imp")
}

