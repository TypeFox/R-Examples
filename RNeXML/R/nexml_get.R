#' Get the desired element from the nexml object
#' 
#' Get the desired element from the nexml object
#' @aliases nexml_get get_item
#' @param nexml a nexml object (from read_nexml)
#' @param element the kind of object desired, see details.  
#' @param ... additional arguments, if applicable to certain elements 
#' @details
#'  
#' \itemize{
#'  \item{"tree"}{ an ape::phylo tree, if only one tree is represented.  Otherwise returns a list of lists of multiphylo trees.  To consistently recieve the list of lists format (preserving the heriarchical nature of the nexml), use \code{trees} instead.}
#'  \item{"trees"}{ returns a list of lists of multiphylo trees, even if all trees are in the same `trees` node (and hence the outer list will be of length 1) or if there is only a single tree (and hence the inner list will also be of length 1.  This guarentees a consistent return type regardless of the number of trees present in the nexml file, and also preserves any heirarchy/grouping of trees.  }
#'  \item{"flat_trees"}{ a multiPhylo object (list of ape::phylo objects) Note that this method collapses any heirachical structure that may have been present as multiple `trees` nodes in the original nexml (though such a feature is rarely used).  To preserve that structure, use `trees` instead.}
#'  \item{"metadata"}{Get metadata from the specified level (default is top/nexml level) }
#'  \item{"otu"}{ returns a named character vector containing all available metadata.  names indicate \code{property} (or \code{rel} in the case of links/resourceMeta), while values indicate the \code{content} (or \code{href} for links). }
#'  \item{"taxa"}{ alias for otu }
#' }
#' For a slightly cleaner interface, each of these elements is also defined as an S4 method
#' for a nexml object.  So in place of `get_item(nexml, "tree")`, one could use `get_tree(nexml)`,
#' and so forth for each element type.  
#' @return return type depends on the element requested.  See details.  
#' @export
#' @seealso \code{\link{get_trees}}
#' @include classes.R
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' nexml_get(nex, "trees")
#' nexml_get(nex, "characters_list")
nexml_get <- function(nexml, 
                     element = c("trees", 
                                 "trees_list", 
                                 "flat_trees", 
                                 "metadata", 
                                 "otu",
                                 "taxa",
                                 "characters", 
                                 "characters_list",
                                 "namespaces"), 
                     ...){
  element <- match.arg(element)

  switch(element,
         trees = get_trees(nexml), # will warn if more than one tree is available
         trees_list = get_trees_list(nexml),
         flat_trees = get_flat_trees(nexml),
         metadata = get_metadata(nexml, ...),
         otu = get_taxa(nexml),
         taxa = get_taxa(nexml),
         characters = get_characters(nexml),
         characters_list = get_characters_list(nexml),
         namespaces = get_namespaces(nexml))
}

get_item <- nexml_get


