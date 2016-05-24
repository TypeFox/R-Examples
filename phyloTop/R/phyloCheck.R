#' Check tree format (internal)
#' 
#' Coerces the tree into a rooted, binary tree of class \code{phylo}. Note that this function used to require trees to be of class \code{phylo4} but the package now uses class \code{phylo} throughout.
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. For most functions in this package the tree should be binary and rooted, hence this function is called to check.
#' If necessary the tree will be coerced into a binary rooted tree using multi2di, if possible.
#' 
#' @return A binary, rooted tree of class \code{phylo}, if possible.
#' 
#' @import ape
#' @importFrom methods as
#' 
#' @keywords internal
#' 
#' @export
phyloCheck <- function(tree) {
  # check input class
  if (class(tree)=="phylo4") {tree <- as(tree, "phylo")}
  else if (class(tree) != "phylo") {stop("Input of class phylo or phylo4 expected.")}
  
  # check if the tree is binary and rooted; if not, do multi2di with a warning
  if(!is.binary.tree(tree)||!is.rooted(tree)) {
    warning("A binary, rooted tree is expected. Applying multi2di to the supplied tree.") 
    tree <- multi2di(tree, random=FALSE)
    # check if it worked (sometimes it can't manage it!)
    if(!is.binary.tree(tree)||!is.rooted(tree)) {stop("Unable to coerce the tree to be binary and rooted.")}
  }
  return(tree)
}