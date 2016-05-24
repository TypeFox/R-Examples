plotccda.cluster <-
function(x){
  par(mfrow=c(1,1))
  plot(x$cluster, main="Cluster Dendrogram according to means", sub="(using squared Euclidean distances and Ward's method)", xlab="groupings", hang=-1 ,font.main=3,labels=x$nameslist)
  x$cluster
}
