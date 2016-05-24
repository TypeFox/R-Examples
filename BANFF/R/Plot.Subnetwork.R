######Sub-network Plot
Plot.Subnetwork=function(net,trace)
{
  ztall=sapply(1:ncol(trace),function(kk) return(mean(trace[,kk])))
  eids=which(ztall>0.5)
  g <- graph.adjacency(net,mode="undirected")
  
  plot(g,layout=layout.kamada.kawai,mark.groups=eids,vertex.size=1,vertex.label=NA)
  
}
