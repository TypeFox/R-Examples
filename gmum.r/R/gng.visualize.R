
.gng.plot2d.errors<-function(gngServer, vertex.color, layout, vertex.size=3){
  ig <- convertToIGraph(gngServer)
  
  if(length(V(ig))==0){
    return()
  }
  
  if(vertex.color == 'label'){
    vertex.color = c(1:length(V(ig)))
    max_col = 0
    for(label in V(ig)$label)
      max_col = max(max_col, round(label))
    cols = rainbow(max_col+1)
    vertex.color = cols[as.double(lapply(V(ig)$label, round))]
  }
  
  if(vertex.color == 'component'){
    vertex.color <- predictComponent(gngServer, )
  }
  
  .visualizeIGraph2dWithErrors(ig, vertex.color, layout, gngServer, vertex.size=3)
}

.gng.plot2d<-function(gngServer, vertex.color, layout, vertex.size=3){
  ig <- convertToIGraph(gngServer)
  
  if(length(V(ig))==0){
    return()
  }
  
  if(vertex.color == 'label'){
    vertex.color = c(1:length(V(ig)))
    max_col = 0
    for(label in V(ig)$data.label)
      max_col = max(max_col, round(label))
    cols = rainbow(max_col+1)
    vertex.color = cols[as.double(lapply(V(ig)$data.label, round))]
  }
  
  .visualizeIGraph2d(ig, vertex.color, layout, vertex.size=vertex.size)
}

# Visualize igraph using igraph plot
# It will layout graph using v0 and v1 coordinates
# @note It is quite slow, works for graphs < 2000 nodes, and for graphs <400 when using layout
.visualizeIGraph2d<-function(g, vertex.color, layout, vertex.size=3){
  L<-layout(g)
  if(vertex.color == 'cluster'){   
    communities <- infomap.community(g)
    communities
    col <- rainbow(length(communities))
    vertex.color <- col[membership(communities)]
  }
  else if(vertex.color == 'fast_cluster'){
    l = fastgreedy.community(g)#as.undirected(g))
    col <- rainbow(length(l))
    print(membership(l))
    vertex.color <- col[membership(l)]
  }
  else if(vertex.color == 'none'){
    vertex.color = NA
  }else{
    # Passed something else as vector
  }
  
  plot.igraph(g,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,layout=L)
}

.visualizeIGraph2dWithErrors<-function(ig, vertex.color, layout_2d, gng,vertex.size=3){
  plot.new()
  par(mfrow=c(1,2))
  .visualizeIGraph2d(ig, vertex.color, layout_2d,vertex.size=vertex.size)
  title("Graph visualization")
  errors_raw = gng$getErrorStatistics()
  errors_raw = errors_raw[5:length(errors_raw)]
  errors = errors_raw
  #errors = log((errors_raw)/min(errors_raw+1e-4))
  plot(errors, type="l", lty=2, lwd=2, xlab="Batch", ylab="Mean batch error", frame.plot=F)
  title("Mean error (log)")
}
