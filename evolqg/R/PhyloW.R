#'Calculates ancestral states of some statistic
#'
#'Calculates weighted average of covariances matrices along a phylogeny, returning a withing-group covariance matrice for each node.
#'
#'@param tree phylogenetic tree
#'@param tip.data list of tip nodes covariance matrices
#'@param tip.sample.size vector of tip nodes sample sizes
#'@return list with calculated within-group matrices, using labels or numbers from tree
#'@export
#'@importFrom ape reorder.phylo
#'@import plyr
#'@examples
#'library(ape)
#'data(dentus)
#'data(dentus.tree)
#'tree <- dentus.tree
#'mat.list <- dlply(dentus, 'species', function(x) cov(x[,1:4]))
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'PhyloW(tree, mat.list, sample.sizes)

PhyloW <- function(tree, tip.data, tip.sample.size = NULL){
  if(is.null(tree$node.label)){
    node.names <- tree$tip.label
  } else{
    node.names <- c(tree$tip.label, tree$node.label)
  }
  if(is.null(tip.sample.size))
    tip.sample.size <- rep(1, length(tip.data))
  names(tip.sample.size) <- names(tip.data)
  tip.sample.size <- as.list(tip.sample.size)
  ancestral.stats <- tip.data
  if(!all(tree$tip.label %in% names(tip.data))) stop("All tip labels must be in stat list.")
  node.order <- unique(reorder(tree, "postorder")$edge[,1])
  for (node in node.order){
      if(is.na(node.names[node])){
        node.names[node] <- as.character(node)
      }
      descendants.list <- node.names[tree$edge[which(tree$edge[,1]==node),2]]
      ponderados <- Map('*', ancestral.stats[descendants.list], tip.sample.size[descendants.list])
      tip.sample.size[[node.names[node]]] <- Reduce("+", tip.sample.size[descendants.list])
      ancestral.stats[[node.names[node]]] <- Reduce("+", ponderados)/
                                              tip.sample.size[[node.names[node]]]
  }
  return(ancestral.stats)
}
