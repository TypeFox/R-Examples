#' Converts Extended Newick Format to Traditional Newick Format
#' 
#' Takes a phylogenetic tree in extended Newick format and converts it to 
#' traditional Newick format that can be directly manipulated by packages such 
#' as \code{ape} and \code{phytools}. 
#' @param eNewick    phylogenetic tree in extended Newick format.
#' @return           phylogenetic tree in traditional Newick format.
#' @export
#' @examples 
#'  ### Converts the phylogenetic tree into traditional Newick format. 
#'  tree <- "((A:0.1,(B:0.3,C:0.2):0.2[60]):0.4[100],(E:0.12,F:0.09):0.4[100]);"
#'  new.tree <- convert.eNewick(tree)
#'  new.tree

convert.eNewick <- function(eNewick) {
  # eNewick = (...):branch_length[node_support],...
  # Newick =  (...)node_support:branch_length,...
  ## Test if eNewick argument is a string.
  if (! is.character(eNewick) || is.list(eNewick) || length(eNewick) > 1) {
    err <- simpleError(call = "eNewick", message = "ERROR: 'eNewick' is not a string!")
    stop(err)
  }
  
  # Split text by ')'; symbolizes a node boundary.
  tree.split <- strsplit(eNewick, ")")
  tree.new <- list()
  # For each split 'node'.
  for (x in 1:length(tree.split[[1]])) {
    # Takes only the support value from eNewick ":x.xxx[<HERE>]..." 
    # If gsub() does not find a match, node.support == tree.split[[1]][x]
    node.support <- gsub(".*\\[(.*)\\].*", "\\1", tree.split[[1]][x])
    # Check that node.support != tree.split[[1]][x] and thus has found the pattern.
    if (! tree.split[[1]][x] %in% node.support) {
      # Takes everything but '<[node_support]>'.
      node.new <- gsub("(.*)\\[.*\\](.*)", "\\1\\2", tree.split[[1]][x])
      # Joins the 'node_support:node.new'.
      tree.new <- c(tree.new, paste(node.support, node.new, sep=""))
    } else { # If no node_support value found.
      tree.new <- c(tree.new, tree.split[[1]][x])
    }
  }
  # Paste converted parts back together.
  Newick <- paste(tree.new, collapse=")")
  return (Newick)
}