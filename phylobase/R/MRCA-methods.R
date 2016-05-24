
##' Most Recent Common Ancestor (MRCA) of 2 or more nodes.
##'
##' Given some nodes (i.e., tips and/or internal), this function
##' returns the node corresponding to the most recent common ancestor.
##'
##' If \code{phy} is a \code{phylo4} or \code{phylo4d} object, the
##' nodes can contain both numeric or character values that will be
##' used by \code{getNode} to retrieve the correct node. However, if
##' \code{phy} is a \code{phylo} object, the nodes must be a numeric
##' vector.
##'
##' With \code{phylo4} and \code{phylo4d} objects, if a single node is
##' provided, it will be returned.
##' 
##' @title MRCA
##' @param phy a phylogenetic tree in phylo4, phylo4d or phylo format.
##' @param ... a vector of nodes
##' @return the node corresponding to the most recent common ancestor
##' @export
##' @include phylo4d-methods.R getNode-methods.R
##' @include oldclasses-class.R
##' @rdname MRCA
##' @examples
##'   data(geospiza)
##'   MRCA(geospiza, 1, 5)
##'   MRCA(geospiza, "fortis", 11)
##'   MRCA(geospiza, 2, 4, "fusca", 3)
##'   geo <- as(geospiza, "phylo")
##'   MRCA(geo, c(1,5))
setGeneric("MRCA", function(phy, ...) {
    standardGeneric("MRCA")
})

##' @rdname MRCA
##' @aliases MRCA,phylo4-method
setMethod("MRCA", signature(phy = "phylo4"), function(phy, ...) {
    nodes <- list(...)
    ## if length==1 and first element is a vector,
    ##   use it as the list
    if (length(nodes)==1 && length(nodes[[1]])>1) {
        nodes <- as.list(nodes[[1]])
    }

    lNodes <- sapply(nodes, function(nd) {
        getNode(x=phy, node=nd, missing="fail")
    })
    
    ## Correct behavior when the root is part of the nodes
    uniqueNodes <- unique(lNodes)
    root <- nodeId(phy, "root")
    if(root %in% uniqueNodes) {
        res <- getNode(phy, root)
        return(res)
    }
    
    ## Correct behavior in case of MRCA of identical taxa
    if(length(uniqueNodes) == 1) {
        res <- uniqueNodes[[1]]
        return(res)
    }
    else {
        ancests <- lapply(nodes, ancestors, phy=phy, type="ALL")
        res <- getNode(phy, max(Reduce(intersect, ancests)))
        return(res)
    }
})

##' @rdname MRCA
##' @aliases MRCA,phylo-method
setMethod("MRCA", signature(phy = "phylo"), function(phy, ...) {
    ape::getMRCA(phy, ...)
})
