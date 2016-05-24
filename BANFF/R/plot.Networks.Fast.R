########print the class "Networks.Fast"
plot.Networks.Fast=function(x, ...){
  ztall=sapply(1:ncol(x$trace),function(kk) return(mean(x$trace[,kk])))
  eids=which(ztall>0.5)
  g <- x$graph
  plot(g,layout=layout.kamada.kawai,mark.groups=eids,vertex.size=1,vertex.label=NA,mark.shape=0)
}
