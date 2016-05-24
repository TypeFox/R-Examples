# force directed graph on all pairs association matrix
fdg <- function(dataSet, dataName=NULL, method="A", cutoff=0.1, dim=2){
  
  d <- dataSet[complete.cases(dataSet),]
  d <- d[,sapply(d,is.numeric)]
  
  # set up adjacency matrix
  if (method == "A") {
    adj <- tap(d)
    adj <- adj + diag(nrow(adj))
    
  } else {
    adj <- cor(d)^2
    adj <- adj - diag(nrow(adj))
  }
  
  # apply cutoff
  adj[adj<cutoff] <- 0
  
  # generate graph
  gr <- graph.adjacency(as.matrix(adj), weighted=TRUE,  mode="upper")
  
  # color the largest cliques green
  gr <- set.vertex.attribute(gr, "color", index=V(gr), "LightBlue")
  gr <- set.vertex.attribute(gr, "color", index=unlist(largest.cliques(gr)), "LightGreen")
  
  if (dim == 2) {
    gr <- set.vertex.attribute(gr, "label", index=V(gr), V(gr)$name)
    gr <- set.edge.attribute(gr, "label", index=E(gr), round(100*E(gr)$weight))
    gr <- set.edge.attribute(gr, "width", index=E(gr), 10 * E(gr)$weight)
  } else {
    gr <- set.vertex.attribute(gr, "label", index=V(gr), paste("   ",V(gr)$name,sep=" "))
    gr <- set.vertex.attribute(gr, "label.color", index=V(gr), "white")
    # the following causes subscript out of bounds in rglplot.igraph
    # gr <- set.edge.attribute(gr, "label", index=E(gr), round(100*E(gr)$weight))
    gr <- set.edge.attribute(gr, "width", index=E(gr), 10 * E(gr)$weight)
  }
  
  # layout <- layout.fruchterman.reingold(gr,dim=dim,coolexp=1)
  layout <- layout.kamada.kawai(gr,dim=dim,coolexp=0.99)
  main <- "Force Directed Graph\n"
  if(!( is.null(dataName))) main <- paste(main, "name ~",dataName,",")
  main <- paste(main,"attraction ~",method,",")
  main <- paste(main,"cutoff ~",100*cutoff,"%")
  
  if(dim == 2) {
    plot.igraph(gr, layout=layout, main=paste("new",main))
  } else {
    rglplot(gr, layout=layout)
  }
}
