ReorderData <- function(tree, data, taxa.names="row names"){
  new.data <- data
  if(is.vector(data)){
    for(i in 1:length(tree$tip.label)){
      new.data[i] <- data[which(names(data) == tree$tip.label[i])]
      names(new.data)[i] <- names(data)[which(names(data) == tree$tip.label[i])]
    }
  }
  if(is.data.frame(data) || is.matrix(data)){
    if(taxa.names == "row names"){
      row.names(new.data) <- 1:length(tree$tip.label)
      for(i in 1:length(tree$tip.label)){
        new.data[i,] <- data[which(row.names(data) == tree$tip.label[i]),]
        row.names(new.data)[i] <- row.names(data)[which(row.names(data) == tree$tip.label[i])]
      }
    }
    if(is.numeric(taxa.names)){
      for(i in 1:length(tree$tip.label)){
        new.data[i,] <- data[which(data[,taxa.names] == tree$tip.label[i]),]
      }
    }
  }
  return(new.data)  
}
