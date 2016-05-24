
#########################################################
### Ordering
#########################################################

##' reordering trees within phylobase objects
##'
##' Methods for reordering trees into various traversal orders
##'
##' The \code{reorder} method takes a \code{phylo4} or \code{phylo4d}
##' tree and orders the edge matrix (i.e. \code{edges(x)}) in the
##' requested traversal order. Currently only two orderings are
##' permitted, and both require rooted trees. In \code{postorder}, a
##' node's descendants come before that node, thus the root, which is
##' ancestral to all nodes, comes last.  In \code{preorder}, a node is
##' visited before its descendants, thus the root comes first.
##'
##' @name reorder-methods
##' @docType methods
##' @param x a \code{phylo4} or \code{phylo4d} object
##' @param order The desired traversal order; currently only
##' \dQuote{preorder} and \dQuote{postorder} are allowed for
##' \code{phylo4} and \code{phylo4d} objects.
##' @param \dots additional optional elements (not in use)
##' @return A \code{phylo4} or \code{phylo4d} object with the edge,
##' label, length and data slots ordered as \code{order}, which is
##' itself recorded in the order slot.
##' @note The \code{preorder} parameter corresponds to
##' \code{cladewise} in the \code{ape} package, and \code{postorder}
##' corresponds (almost) to \code{pruningwise}.
##'
##' @author Peter Cowan, Jim Regetz
##' @seealso \code{\link[ape]{reorder.phylo}} in the \code{ape} package.
##' \code{\link{ancestors}} \code{\link{ancestor}} \code{\link{siblings}}
##' \code{\link{children}} \code{\link{descendants}}
##' @keywords methods
##' @include phylo4-class.R
##' @include phylo4-methods.R
##' @exportMethod reorder
##' @aliases reorder
##' @examples
##' phy <- phylo4(ape::rtree(5))
##' edges(reorder(phy, "preorder"))
##' edges(reorder(phy, "postorder"))
setGeneric("reorder")

##' @rdname reorder-methods
##' @aliases reorder,phylo4-method
setMethod("reorder", signature(x="phylo4"),
 function(x, order=c("preorder", "postorder")) {
    ## call orderIndex and use that index to order edges, labels and lengths
    order   <- match.arg(order)
    index   <- orderIndex(x, order)
    x@order <- order
    x@edge  <- edges(x)[index, ]
    if(hasEdgeLabels(x)) {
        x@edge.label  <- x@edge.label[index]
    }
    if(hasEdgeLength(x)) {
        x@edge.length <- x@edge.length[index]
    }
    x
})

## non exported function
orderIndex <- function(x, order=c("preorder", "postorder")) {

    order <- match.arg(order)
    if(!isRooted(x)){
        stop("Tree must be rooted to reorder")
    }
    ## get a root node free edge matrix
    edge <- edges(x, drop.root=TRUE)
    ## Sort edges -- ensures that starting order of edge matrix doesn't
    ## affect the order of reordered trees
    edge <- edge[order(edge[, 2]), ]

    # recast order argument as integer to pass to C
    if(order == 'postorder') {
        iOrder <- 0L
    } else if(order == 'preorder') {
        iOrder <- 1L
    } else {stop(paste("Method for", order, "not implemented"))}

    if (!hasPoly(x) & !hasSingle(x)) {
        # method 1: faster, but only works if all internal nodes have
        # exactly two children (true binary tree)

        # extract nodes, separating descendants into left (first
        # encountered) and right (second encountered) for each ancestor
        isFirst <- !duplicated(edge[, 1])
        ancestor <- as.integer(edge[isFirst, 1])
        left <- as.integer(edge[isFirst, 2])
        right <- as.integer(edge[!isFirst, 2])[match(ancestor,
            edge[!isFirst, 1])]
        descendantNew <- rep(0L, nEdges(x))
        root <- as.integer(rootNode(x))
        nEdge <- as.integer(length(ancestor))

        descendantReord <- .C("reorderBinary", descendantNew, root,
            ancestor, left, right, nEdge, iOrder)[[1]]

    } else {
        # method 2: not as fast, but robust to singletons and polytomies

        # extract ancestors and descendants
        ancestor <- as.integer(edge[,1])
        descendant <- as.integer(edge[,2])
        descendantNew <- rep(0L, nEdges(x))
        root <- as.integer(rootNode(x))
        nEdge <- as.integer(nrow(edge))

        descendantReord <- .C("reorderRobust", descendantNew, root,
            ancestor, descendant, nEdge, iOrder)[[1]]

    }

    ## Original pure R implementation of the above:
    #### recursive functions are placed first and calls to those functions below
    ##postOrder <- function(node) {
    ##    ## this function returns a vector of nodes in the post order traversal
    ##    ## get the descendants
    ##    traversal <- NULL
    ##    ## edge -- defined above, outside this function
    ##    ## extensive testing found this loop to be faster than apply() etc.
    ##    for(i in edge[edge[, 1] == node, 2]) {
    ##        traversal <- c(traversal, postOrder(i))
    ##    }
    ##    c(traversal, node)
    ##}
    ##preOrder  <- function(node) {
    ##    ## see expanded code in comments of postOrder()
    ##    ## only difference here is that we record current node, then descendants
    ##    traversal <- NULL
    ##    for(i in edge[edge[, 1] == node, 2]) {
    ##        traversal <- c(traversal, preOrder(i))
    ##    }
    ##    c(node, traversal)
    ##}
    ##if(order == 'postorder') {
    ##    descendantReord <- postOrder(rootNode(x))
    ##} else if(order == 'preorder') {
    ##    descendantReord <- preOrder(rootNode(x))
    ##} else {stop(paste("Method for", order, "not implemented"))}

    ## match the new node order to the old order to get an index
    index <- match(descendantReord, edges(x)[, 2])

}
