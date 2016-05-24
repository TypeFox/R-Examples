distmap <-
function(data, what=c("subjects", "variables"), col, ...){
#function(data, what=c("subjects", "variables"), type=list(), col, ...){
  # !! to be done: allow also asymmetric binary variables

  what <- match.arg(what)
  
  if(what == "subjects"){
    D <- dist.subjects(data)
    #D <- dist.subjects(data, type=type)

    # similarity matrix
    S <- 1 - D^2
    
    # corresponding dendrogram
    dendro <- as.dendrogram(hclust(D))
  }
  else if(what == "variables"){
    S <- similarity.variables(data)
    dendro <- as.dendrogram(hclust(as.dist(sqrt(1 - S))))
  }
  
  # graphical parameters
  if(missing(col))
    col <- marray::maPalette(low="#F7FBFF", mid="#6BAED6", high="#08306B", 30)

  # distogram
  gplots::heatmap.2(as.matrix(S), Rowv=dendro, Colv=dendro, trace="none", density.info="none", col=col, ...)
}
