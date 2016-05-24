FuzzyMatch <- function(tree, data, max.dist){
  tree.names <- tree$tip.label    # lets get the names on the tree
  close.taxa <- data.frame()
  counter <- 1
  data.names <- unique(data)
  for(i in 1:length(data.names)){
    name.dist <- min(adist(data.names[i], as.character(tree.names)))
    if( name.dist <= max.dist & name.dist > 0){
      fuq <- which.min(adist(data.names[i], as.character(tree.names)))
      close.taxa[counter, 1] <- data.names[i]
      close.taxa[counter, 2] <- tree.names[fuq]  
      close.taxa[counter, 3] <- name.dist
      counter <- counter + 1
    }
  }
  if(counter == 1){
    cat("Found", counter - 1, "names that were close but imperfect matches\n")
  }
  if(counter > 1){
    cat("Found", counter - 1, "names that were close but imperfect matches\n")
    colnames(close.taxa) <- c("name.in.data", "name.in.tree", "differences")
    return(close.taxa)
  }
}

