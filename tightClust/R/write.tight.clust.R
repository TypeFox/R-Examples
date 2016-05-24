write.tight.clust <-
function(x, ...) {
    data<-x$data
    cluster<-matrix(x$cluster,ncol=1)
    colnames(cluster)<-"Cluster"
    data<-cbind(data,cluster)
    order.genes<-unlist(lapply(c(1:max(cluster),-1),function(i) which(cluster==i)))
    data<-data[order.genes,]
    write.table(data, ...)
}
