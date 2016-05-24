#' Split topology
#'
#' For each node at a given distance from the root, this function finds the size of its induced subclade, i.e. its number of tip descendants.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Michael Boyd \email{mboyd855@gmail.com}
#'
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param dist integer distance of nodes of interest from the root.
#' 
#' @return A vector of integers, each corresponding to the clade size of a node at the given distance from the root. The clade sizes are given in ascending order and each is labelled by its node name or number.
#' This vector can be considered as a partition of the tips or the "split topology" of the tree at a given depth.
#' 
#' @import ape
#'
#' @examples
#' ## Find the split topology of a random tree with 20 tips, at a distance 2 from the root:
#' tree <- rtree(20)
#' plot(tree)
#' splitTop(tree,2)
#' 
#'
#' @export
splitTop <- function(tree,dist) {
  tree <- phyloCheck(tree)
  depths <- getDepths(tree)
  allDepths <- c(depths$tipDepths,depths$nodeDepths)
  if (dist > max(allDepths)) {stop(paste0("For this tree, 'dist' must be <=",max(allDepths),", the maximum node depth in the tree."))}
  
  # IDs of nodes at given distance:
  nodes <- which(allDepths==dist)
  
  # find all clade sizes:
  configs <- nConfig(tree)$cladeSizes
  # extract clade sizes of nodes at given distance:
  splits <- configs[nodes]
  
  return(sort(splits))
}
