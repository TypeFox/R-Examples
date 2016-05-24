
##' nodeId methods
##'
##' These functions gives the node (\code{nodeId}) or edge
##' (\code{edgeId}) identity.
##'
##' \code{nodeId} returns the node in ascending order, and
##' \code{edgeId} in the same order as the edges are stored in the
##' edge matrix.
##' 
##' @param x a \code{phylo4} or \code{phylo4d} object.
##' @param type a character vector indicating which subset of the
##' nodes or edges you are interested in.
##' @return \describe{
##'   \item{nodeId}{an integer vector  indicating node numbers}
##'   \item{edgeId}{a character vector indicating the edge identity}
##' }
##' @export
##' @docType methods
##' @rdname nodeId-methods
##' @include phylo4-class.R phylo4-methods.R phylo4-accessors.R root-methods.R
##' @examples
##'   data(geospiza)
##'   identical(nodeId(geospiza, "tip"), 1:nTips(geospiza))
##'   nodeId(geospiza, "internal")
##'   edgeId(geospiza, "internal")
##'   nodeId(geospiza, "root")
setGeneric("nodeId", function(x, type=c("all", "tip", "internal",
    "root")) {
    standardGeneric("nodeId")
})

##' @rdname nodeId-methods
##' @aliases nodeId,phylo4-method
setMethod("nodeId", signature(x="phylo4"),
 function(x, type=c("all",
    "tip","internal","root")) {

     type <- match.arg(type)
     E <- edges(x)

     ## Note: this implementation will still work even if tips are not
     ## 1:nTips and nodes are not (nTips+1):nNodes
     nid <- switch(type,
         ## all nodes appear at least once in the edge matrix
         ## twice slower: all = unique(as.vector(E)[as.vector(E) != 0]),
         ## but maybe should be used if tree is not "normal"
         all = {
             if (isRooted(x)) {
                 res <- getAllNodesFast(x@edge)[-1]
             }
             else {
                 res <- getAllNodesFast(x@edge)
             }
             res
         },
         ## tips are nodes that do not appear in the ancestor column
         ## three times slower: setdiff(E[, 2], E[, 1]),
         tip = tipsFast(x@edge[,1]),
         ## internals are nodes that *do* appear in the ancestor column
         ## about 0.5 faster than: setdiff(getAllNodesFast(x@edge, isRooted(x)), tipsFast(x@edge[,1])),
         internal = unique(E[E[, 1] != 0, 1]), 
         ## roots are nodes that have NA as ancestor
         root = if (!isRooted(x)) return(NA) else unname(E[E[, 1] == 0, 2]))

     return(sort(nid))

})

#### ----- edgeId

##' @rdname nodeId-methods
##' @aliases edgeId
##' @export
setGeneric("edgeId", function(x, type=c("all", "tip", "internal",
    "root")) {
    standardGeneric("edgeId")
})

##' @rdname nodeId-methods
##' @aliases edgeId,phylo4-method
setMethod("edgeId", signature(x="phylo4"),
 function(x, type=c("all", "tip",
    "internal", "root")) {
    type <- match.arg(type)
    edge <- edges(x)
    if (type=="tip") {
        isTip <- !(edge[, 2] %in% edge[, 1])
        edge <- edge[isTip, , drop=FALSE]
    } else if (type=="internal") {
        isInt <- (edge[, 2] %in% edge[, 1])
        edge <- edge[isInt, , drop=FALSE]
    } else if (type=="root") {
        isRoot <- edge[, 1] == 0
        edge <- edge[isRoot, , drop=FALSE]
    } # else just use complete edge matrix if type is "all"
    id <- paste(edge[, 1], edge[, 2], sep="-")
    return(id)
})

