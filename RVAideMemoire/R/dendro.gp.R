dendro.gp <- function(dend) {
  d <- rev(diff(dend$height))
  k <- which.max(d)+1
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  par(mfrow=c(1,2))
  plot(dend,hang=-1)
  rect.hclust(dend,k=k)
  xlab <- paste("Number of groups\n(best value from dendrogram: ",k,")",sep="")
  g <- barplot(d,xlab=xlab,ylab="Difference in agglomeration index with the next level",
    names.arg=2:length(dend$height))
  demi <- (g[2]-g[1])/2
  rect(g[k-1]-demi,0,g[k-1]+demi,d[k-1],border="red")  
}
