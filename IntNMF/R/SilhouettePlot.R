SilhouettePlot <-
function(fit,cluster.col=NULL){
distance <- 1-fit$consensus
if (is.null(cluster.col)) cluster.col <- c("gray",2:length(unique(fit$clusters)))
if (length(unique(fit$clusters))!=length(cluster.col)) stop("Specified number of colors didn't match with number of clusters")
plot(silhouette(fit$clusters,dmatrix=distance),col=cluster.col,main="Silhouette plot")
print(summary(silhouette(fit$clusters,dmatrix=distance)))
cat("------------------------------------------------\n")
cat("Mean silhouette width = ",mean((silhouette(fit$clusters,dmatrix=distance))[,3]),"\n")
}
