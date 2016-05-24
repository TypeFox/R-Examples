

.shortestPathInt <- function(phy, node1, node2){
    ## some checks
    ## if (is.character(checkval <- checkPhylo4(x))) stop(checkval) # no need
    t1 <- getNode(phy, node1)
    t2 <- getNode(phy, node2)
    if(any(is.na(c(t1,t2)))) stop("wrong node specified")
    if(t1==t2) return(NULL)

    ## main computations
    comAnc <- MRCA(phy, t1, t2) # common ancestor
    desComAnc <- descendants(phy, comAnc, type="all")
    ancT1 <- ancestors(phy, t1, type="all")
    path1 <- intersect(desComAnc, ancT1) # path: common anc -> t1

    ancT2 <- ancestors(phy, t2, type="all")
    path2 <- intersect(desComAnc, ancT2) # path: common anc -> t2

    res <- union(path1, path2) # union of the path
    ## add the common ancestor if it differs from t1 or t2
    if(!comAnc %in% c(t1,t2)){
        res <- c(comAnc,res)
    }

    res <- getNode(phy, res)

    return(res)
}

##' Finds the shortest path between two nodes in a tree
##'
##' Given two nodes (i.e, tips or internal nodes), this function
##' returns the shortest path between them (excluding \code{node1} and
##' \code{node2} as a vector of nodes.
##' @title shortestPath-methods
##' @param x a tree in the phylo4, phylo4d or phylo format
##' @param node1 a numeric or character (passed to \code{getNode})
##' indicating the beginning from which the path should be calculated.
##' @param node2 a numeric or character (passed to \code{getNode})
##' indicating the end of the path.
##' @return a vector of nodes indcating the shortest path between 2 nodes
##' @seealso getNode
##' @rdname shortestPath-methods
##' @docType methods
##' @include MRCA-methods.R
##' @export
setGeneric("shortestPath", function(x, node1, node2) {
    standardGeneric("shortestPath")
})

##' @name shortestPath-phylo4
##' @rdname shortestPath-methods
##' @aliases shortestPath,phylo4-method
setMethod("shortestPath", signature(x="phylo4", node1="ANY", node2="ANY"),
          function(x, node1, node2) {
              .shortestPathInt(phy=x, node1=node1, node2=node2)
          })

##' @name shortestPath-phylo
##' @rdname shortestPath-methods
##' @aliases shortestPath,phylo-method
setMethod("shortestPath", signature(x="phylo", node1="ANY", node2="ANY"),
          function(x, node1, node2) {
              phy <- as(x, "phylo4")
              .shortestPathInt(phy=phy, node1=node1, node2=node2)
          })
