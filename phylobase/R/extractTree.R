## extract the phylo4 part of phylo4d; relies on implicit coerce method

##' Get tree from tree+data object
##' 
##' Extracts a \code{phylo4} tree object from a \code{phylo4d}
##' tree+data object.
##' 
##' \code{extractTree} extracts just the phylogeny from a tree+data
##' object. The phylogeny contains the topology (how the nodes are
##' linked together), the branch lengths (if any), and any tip and/or
##' node labels. This may be useful for extracting a tree from a
##' \code{phylo4d} object, and associating with another phenotypic
##' dataset, or to convert the tree to another format.
##' 
##' @param from a \code{phylo4d} object, containing a phylogenetic
##' tree plus associated phenotypic data. Created by the
##' \code{phylo4d()} function.
##' @author Ben Bolker
##' @seealso \code{\link{phylo4-methods}},
##' \code{\link{phylo4d-methods}}, \code{\link{coerce-methods}} for
##' translation functions.
##' @keywords methods
##' @export
##' @include setAs-methods.R
##' @examples
##' tree.phylo <- ape::read.tree(text = "((a,b),c);")
##' tree <- as(tree.phylo, "phylo4")
##' plot(tree)
##' tip.data <- data.frame(size = c(1, 2, 3), row.names = c("a", "b", "c"))
##' (treedata <- phylo4d(tree, tip.data))
##' plot(treedata)
##' (tree1 <- extractTree(treedata))
##' plot(tree1)
##' 
extractTree <- function(from) {
    as(from, "phylo4")
}
