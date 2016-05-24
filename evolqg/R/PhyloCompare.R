#'Compares sister groups
#'
#'Calculates the comparison of some statistic between sister groups along a phylogeny
#'
#'@param tree phylogenetic tree
#'@param node.data list of node data
#'@param ComparisonFunc comparison function, default is PCAsimilarity
#'@param ... Aditional arguments passed to ComparisonFunc
#'@return list with a data.frame of calculated comparisons for each node, using labels or numbers from tree; and a list of comparisons for plotting using phytools (see examples)
#'@note Phylogeny must be fully resolved
#'@author Diogo Melo
#'@export
#'@importFrom ape reorder.phylo
#'@import plyr
#'@examples
#'library(ape)
#'data(bird.orders)
#'tree <- bird.orders
#'mat.list <- RandomMatrix(5, length(tree$tip.label))
#'names(mat.list) <- tree$tip.label
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'phylo.state <- PhyloW(tree, mat.list, sample.sizes)
#'
#'phylo.comparisons <- PhyloCompare(tree, phylo.state)
#'
#'# plotting results on a phylogeny:
#'library(phytools)
#'plotBranchbyTrait(tree, phylo.comparisons[[2]])
PhyloCompare <- function(tree, node.data, ComparisonFunc = PCAsimilarity, ...){
  if(is.null(tree$node.label)){
    node.names <- tree$tip.label
  } else{
    node.names <- c(tree$tip.label, tree$node.label)
  }
  node.order <- unique(reorder(tree, "postorder")$edge[,1])
  phylo.comparisons <- list()
  tree.comparisons <- list()
  for (node in node.order){
    if(is.na(node.names[node])){
      node.names[node] <- as.character(node)
    }
    current.nodes <- tree$edge[which(tree$edge[,1]==node),2]
    phylo.comparisons[[node.names[node]]] <- ComparisonFunc(node.data[[current.nodes[1]]], node.data[[current.nodes[2]]], ...)
    tree.comparisons[[node.names[current.nodes[1]]]] <- phylo.comparisons[[node.names[node]]][1]
    tree.comparisons[[node.names[current.nodes[2]]]] <- phylo.comparisons[[node.names[node]]][1]
    names(tree.comparisons)[which(names(tree.comparisons) == node.names[current.nodes[1]])] <- paste(node, current.nodes[1], sep = ',')
    names(tree.comparisons)[which(names(tree.comparisons) == node.names[current.nodes[2]])] <- paste(node, current.nodes[2], sep = ',')
  }
  tree.order = aaply(tree$edge, 1, function(x) paste(x, collapse = ','))
  phylo.comparisons.df = ldply(phylo.comparisons, .id = 'node')
  return(list(df = phylo.comparisons.df, tree = tree.comparisons[tree.order]))
}
