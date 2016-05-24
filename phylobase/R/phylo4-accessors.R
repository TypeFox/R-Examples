
##' Number of tips, nodes and edges found in a tree.
##'
##' Function to return the number of tips, nodes and edges found in a
##' tree in the \code{phylo4} or \code{phylo4d} format.
##' @title nTips, nNodes, nEdges
##' @aliases nTips
##' @param x a \code{phylo4} or \code{phylo4d} object
##' @return a numeric vector indicating the number of tips, nodes or
##' edge respectively.
##' @docType methods
##' @export
##' @include phylo4-class.R phylo4-methods.R
##' @include oldclasses-class.R
##' @rdname nTips-methods
setGeneric("nTips", function(x) {
    standardGeneric("nTips")
})

##' @rdname nTips-methods
##' @aliases nTips,phylo4-method
setMethod("nTips", signature(x="phylo4"), function(x) {
    E <- edges(x)
    if(nrow(E) == 0)
        return(0)
    else {
        ## at this time NAs are not allowed in edge matrix
        ## sum(tabulate(E[, 1]) == 0)
        nTipsFastCpp(E[, 1])
    }        
})

##' @rdname nTips-methods
##' @aliases nTips,phylo-method
setMethod("nTips", signature(x="phylo"),
 function(x) {
     Ntip(x)
})

##' @rdname nTips-methods
##' @aliases nNodes
##' @export
setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})

##' @rdname nTips-methods
##' @aliases nNodes,phylo4-method
setMethod("nNodes", signature(x="phylo4"), function(x) {
    E <- edges(x, drop.root=TRUE)
    if(nrow(E) == 0) {
        return(0)
    } else {
        return(length(unique(E[, 1])))
    }
})

##' @rdname nTips-methods
##' @aliases nEdges
##' @export
setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})

##' @rdname nTips-methods
##' @aliases nEdges,phylo4-method
setMethod("nEdges", signature(x="phylo4"),
 function(x) {
    nrow(x@edge)
})


#########################################################
### Edge accessors
#########################################################

##' Edges accessors
##'
##' Access or modify information about the edges.
##'
##' @param x a \code{phylo4} or \code{phylo4d} object.
##' @param drop.root logical (default FALSE), should the edge
##' connecting the root be included in the edge matrix?
##' @param \dots Optional arguments used by specific methods. (None
##' used at present).
##' @return \describe{
##' \item{\code{edges}}{returns the edge matrix that represent the
##' ancestor-descendant relationships among the nodes of the tree.}
##'
##' \item{\code{edgeOrder}}{returns the order in which the edge matrix
##' is in.}
##'
##' \item{\code{internalEdges}}{returns a logical vector indicating
##' internal edges (edges that connect an internal node to
##' another). This vector is named with the \code{edgeId}}.
##'
##' \item{\code{terminalEdges}}{returns a logical vector indicating
##' terminal edges (edges that connect an internal node to a
##' tip). This vector is named with the \code{edgeId} }}
##' @author Ben Bolker, Francois Michonneau, Thibaut Jombart
##' @seealso reorder, edgeId
##' @examples
##'    data(geospiza)
##'    edges(geospiza)
##'    edgeOrder(geospiza)
##'    geoPost <- reorder(geospiza, "postorder")
##'    edgeOrder(geoPost)
##'    ## with a binary tree this should always be true
##'    identical(!terminalEdges(geospiza), internalEdges(geospiza))
##' @export
##' @docType methods
##' @rdname edges-accessors
##' @include phylo4-methods.R
setGeneric("edges", function(x, ...) {
    standardGeneric("edges")
})

##' @rdname edges-accessors
##' @aliases edges,phylo4-method
setMethod("edges", signature(x="phylo4"),
 function(x, drop.root=FALSE) {
     e <- x@edge
     if (drop.root) e <- e[e[, 1] != 0, ]
     e
})

##### -------- edgeOrder

##' @rdname edges-accessors
##' @aliases edgeOrder
##' @export
setGeneric("edgeOrder", function(x, ...) {
  standardGeneric("edgeOrder")
})

##' @rdname edges-accessors
##' @aliases edgeOrder,phylo4-method
setMethod("edgeOrder", signature(x="phylo4"),
 function(x) {
    x@order
})

##### -------- internalEdges

##' @rdname edges-accessors
##' @aliases internalEdges
##' @export
setGeneric("internalEdges", function(x) {
    standardGeneric("internalEdges")
})

##' @rdname edges-accessors
##' @aliases internalEdges,phylo4-method
setMethod("internalEdges", signature(x="phylo4"),
  function(x) {
      res <- edges(x)[, 2] %in% nodeId(x, "internal")
      names(res) <- edgeId(x, "all")
      res
})

##### -------- terminalEdges

##' @rdname edges-accessors
##' @aliases terminalEdges
##' @export
setGeneric("terminalEdges", function(x) {
    standardGeneric("terminalEdges")
})

##' @rdname edges-accessors
##' @aliases terminalEdges,phylo4-method
setMethod("terminalEdges", signature(x="phylo4"),
  function(x) {
      res <- edges(x)[, 2] %in% nodeId(x, "tip")
      names(res) <- edgeId(x, "all")
      res
})
