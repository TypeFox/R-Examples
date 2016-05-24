#' Configuration sizes in a tree
#' 
#' Finds the sizes of configurations in the tree. 
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param maxClade An integer between 1 and the number of tips (the default), specifying the maximum clade size of interest.
#' 
#' @return A list with 2 entries: 
#' \itemize{
#' \item cladeSizes is a vector giving the size of the clade descending at each node. Tips all have the value 1. 
#' Internal nodes have their number of tip descendants. 
#' \item numClades is a vector where numClades[[i]] is the number of clades of size i in the tree. 
#' All clade sizes are calculated, but the output can be restricted using \code{maxClade} to just those of size up to 'maxClade'.
#' }
#' 
#' @import ape
#'  
#' @examples
#' ## Configuration sizes on a random tree with 10 tips:
#' tree <- rtree(10)
#' plot(tree)
#' nodelabels()
#' nConfig(tree)
#' 
#' 
#' @export
nConfig <- function(tree,maxClade=NULL) { 
  tree <- phyloCheck(tree)
  if (is.null(tree$tip.label))
    stop("This tree has no tip labels")
  num.tips=length(tree$tip.label)
  if (is.null(maxClade)) {maxClade <- num.tips}
  else if(maxClade > num.tips) {maxClade <- num.tips} # maxClade greater than number of tips makes no sense and would append unnecessary zeroes to output
  labels <- rep(NA, nrow(tree$edge)+1)
  names(labels)[1:num.tips]=tree$tip.label;
  names(labels)[(num.tips+1): length(labels)]=paste0("node",(num.tips+1): length(labels))
  labels[1:num.tips]=1     # tips are 1 here. 
  NodeIDS= (num.tips + 1) : (2*num.tips -1)
  # fill in the configuration sizes of internal nodes
  while (any(is.na(labels))) { 
    IsReady = NodeIDS[ vapply(NodeIDS,function(x) !any(is.na(labels[tree$edge[which(tree$edge[,1]==x),2]])) & is.na(labels[x])  ,FUN.VALUE=TRUE) ]
    TheseLabels = unlist(sapply(IsReady, function(x) sum(labels[tree$edge[tree$edge[,1]==x,2]])))
    labels[IsReady]=TheseLabels
  }
  numClades=vapply(1:maxClade, function(x) sum(labels==x),FUN.VALUE=1)
  names(numClades) <- 1:maxClade
  return(list(cladeSizes=labels,numClades=numClades))
}