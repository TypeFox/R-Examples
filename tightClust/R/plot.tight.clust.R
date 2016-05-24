plot.tight.clust <-
function(x, standardize.gene=TRUE, order.sample=FALSE, plot.noise=TRUE, ...) {

    mirror.matrix <- function(x) { 
        xx <- as.data.frame(x); 
        xx <- rev(xx); 
        xx <- as.matrix(xx); 
        xx; 
    } 

    rotate180.matrix <- function(x) { 
        xx <- rev(x); 
        dim(xx) <- dim(x); 
        xx; 
    } 

    flip.matrix <- function(x) { 
        mirror.matrix(rotate180.matrix(x)) 
    } 
    
    data<-x$data
    cluster<-x$cluster
    size<-x$size
    if(standardize.gene) data<-t(scale(t(data)))
    l<-1:max(cluster)
    if(plot.noise) l<-c(l,-1)
    order.genes<-unlist(lapply(l,function(i) which(cluster==i)))
    data.plot<-data[order.genes,]
    if(order.sample) data.plot<-data.plot[,hclust(dist(t(data.plot)))$order]
    cuts<-cumsum(size)
    if(!plot.noise) cuts<-cuts[-length(size)]
    nr<-dim(data.plot)[1]
    nc<-dim(data.plot)[2]
    cexCol<-1/log10(nc)
    cexRow<-1/log(nr,4)
    image(x=1:nc,y=1:nr,z=t(flip.matrix(data.plot)),col=colorRampPalette(c("green","black","red"))(64),axes=FALSE,xlab = "", ylab = "", ...)
    axis(4,nr:1,labels=(row.names(data.plot)), las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    for(i in cuts) {
        lines(c(-1,nc+1),rep(nr-i+0.5,2),col="white")
    }
}
