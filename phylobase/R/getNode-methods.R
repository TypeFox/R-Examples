## matching node labels with node numbers ...
## e.g.
## 14 tips, 13 int nodes
## N04 = nodeLabels[4]
##   <-> node 18
## x = n-nTips(phy)
## so:     n = x+nTips(phy)


##' Node and Edge look-up functions
##' 
##' Functions for retrieving node and edge IDs (possibly with corresponding
##' labels) from a phylogenetic tree.
##' 
##' \code{getNode} and \code{getEdge} are primarily intended for looking up the
##' IDs either of nodes themselves or of edges associated with those nodes. Note
##' that they behave quite differently. With \code{getNode}, any input nodes are
##' looked up against tree nodes of the specified type, and those that match are
##' returned as numeric node IDs with node labels (if they exist) as element
##' names. With \code{getEdge}, any input nodes are looked up against edge ends
##' of the specified type, and those that match are returned as character edge
##' IDs with the corresponding node ID as element names.
##' 
##' If \code{missing} is \dQuote{warn} or \dQuote{OK}, \code{NA} is returned for
##' any nodes that are unmatched for the specified type. This can provide a
##' mechanism for filtering a set of nodes or edges.
##' 
##' \code{nodeId} provides similar output to \code{getNode} in the case when no
##' node is supplied, but it is faster and returns an unnamed vector of the
##' numeric IDs of all nodes of the specified node type.  Similarly,
##' \code{edgeId} simply returns an unnamed vector of the character IDs of all
##' edges for which the descendant node is of the specified node type.
##' 
##' @param x a \linkS4class{phylo4} object (or one inheriting from
##' \linkS4class{phylo4}, e.g. a \linkS4class{phylo4d} object)
##' @param node either an integer vector corresponding to node ID numbers, or a
##' character vector corresponding to node labels; if missing, all nodes
##' appropriate to the specified type will be returned by \code{getNode}, and
##' all edges appropriate to the specified type will be returned by
##' \code{getEdge}.
##' @param type (\code{getNode}) specify whether to return nodes matching "all"
##' tree nodes (default), only "tip" nodes, or only "internal" nodes;
##' (\code{nodeId, edgeId}) specify whether to return "all" tree nodes, or only
##' those corresponding to "tip", "internal", or "root" nodes; (\code{getEdge})
##' specify whether to look up edges based on their descendant node
##' ("descendant") or ancestral node ("ancestor")
##' @param missing what to do if some requested node IDs or names are not in the
##' tree: warn, do nothing, or stop with an error
##' @return \item{list("getNode")}{returns a named integer vector of node IDs,
##' in the order of input nodes if provided, otherwise in nodeId order}
##' \item{list("getEdge")}{returns a named character vector of edge IDs, in the
##' order of input nodes if provide, otherwise in nodeId order}
##' \item{list("nodeId")}{returns an unnamed integer vector of node IDs, in
##' ascending order} \item{list("getEdge")}{returns an unnamed character vector
##' of edge IDs, in edge matrix order}
##' @keywords misc
##' @export
##' @rdname getNode-methods
##' @include phylo4-class.R
##' @examples
##' 
##'   data(geospiza)
##'   nodeLabels(geospiza) <- LETTERS[1:nNodes(geospiza)]
##'   plot(as(geospiza, "phylo4"), show.node.label=TRUE)
##'   getNode(geospiza, 18)
##'   getNode(geospiza, "D")
##'   getEdge(geospiza, "D")
##'   getEdge(geospiza, "D", type="ancestor")
##' 
##'   ## match nodes only to tip nodes, flagging invalid cases as NA
##'   getNode(geospiza, c(1, 18, 999), type="tip", missing="OK")
##' 
##'   ## get all edges that descend from internal nodes
##'   getEdge(geospiza, type="ancestor")
##' 
##'   ## identify an edge from its terminal node
##'   getEdge(geospiza, c("olivacea", "B", "fortis"))
##'   getNode(geospiza, c("olivacea", "B", "fortis"))
##'   edges(geospiza)[c(26, 1, 11),]
##' 
##'   ## quickly get all tip node IDs and tip edge IDs
##'   nodeId(geospiza, "tip")
##'   edgeId(geospiza, "tip")
##'
setGeneric("getNode", function(x, node, type=c("all", "tip", "internal"),
                               missing=c("warn", "OK", "fail")) {
    standardGeneric("getNode")
})

##' @rdname getNode-methods
##' @aliases getNode,phylo4-method
setMethod("getNode", signature(x="phylo4", node="ANY"),
          function(x, node, type=c("all", "tip", "internal"),
                   missing=c("warn","OK","fail")) {

    type <- match.arg(type)
    missing <- match.arg(missing)

    ## if missing node arg, get all nodes of specified type
    if (missing(node)) {
        node <- nodeId(x, type)
    }

    if (length(node) == 0) {
      rval <- integer(0)
      names(rval) <- character(0)
      return(rval)
    }
    
    lblTmp <- labels(x, type)
    
    ## match node to tree
    if (is.character(node)) {
        ndTmp <- paste("^\\Q", node, "\\E$", sep="")        
        irval <- lapply(ndTmp, function(ND) {
            grep(ND, lblTmp, perl=TRUE)
        })
        irvalL <- sapply(irval, length)
        irval[irvalL == 0] <- 0
        irval <- unlist(irval)
    } else if (is.numeric(node) && all(floor(node) == node, na.rm=TRUE)) {
        irval <- match(as.character(node), names(lblTmp))
    } else {
        stop("Node must be a vector of class \'integer\' or \'character\'.")
    }

    ## node numbers
    rval <- names(lblTmp)[irval]
    rval[is.na(node)] <- NA # return NA for any NA_character_ inputs, not needed but ensure rval has correct length
    rval <- as.integer(rval)

    ## node labels
    nmNd <- lblTmp[irval]
    names(rval) <- nmNd
    
    ## deal with nodes that don't match
    if (any(is.na(rval))) {
        missnodes <- node[is.na(rval)]
        msg <- paste("Some nodes not found among", type, "nodes in tree:",
            paste(missnodes,collapse=", "))
        if (missing=="fail") {
            stop(msg)
        } else if (missing=="warn") {
            warning(msg)
        }
    }
    return(rval)
})

##' @rdname getNode-methods
##' @aliases getEdge-methods
##' @export
setGeneric("getEdge", function(x, node, type=c("descendant", "ancestor"),
                               missing=c("warn", "OK", "fail")) {
    standardGeneric("getEdge")
})

##' @name getEdge
##' @rdname getNode-methods
##' @aliases getEdge,phylo4-method
setMethod("getEdge", signature(x="phylo4", node="ANY"),
          function(x, node, type=c("descendant", "ancestor"),
                   missing=c("warn", "OK", "fail")) {

    type <- match.arg(type)
    missing <- match.arg(missing)
    if (missing(node)) {
        if (type=="descendant") {
            node <- nodeId(x, "all")
        } else if (type=="ancestor") {
            node <- nodeId(x, "internal")
        }
    }

    node.id <- getNode(x, node, missing="OK")

    nd <- lapply(node.id, function(nid) {
        if (is.na(nid)) {
            res <- NA
        } else {
            res <- switch(type,
                descendant = edgeId(x)[edges(x)[,2] %in% nid],
                ancestor = edgeId(x)[edges(x)[,1] %in% nid])
            ## hack to return NA for tip nodes when type='ancestor'
            if(length(res)==0) res <- NA
            names(res) <- rep(nid, length(res))
        }
        names(res) <- rep(nid, length(res))
        res
    })

    ## warn or stop if necessary
    is.missing <- is.na(nd)
    if (missing!="OK" && any(is.missing)) {
        msg <- paste("Not all nodes are ", type, "s in this tree: ",
            paste(node[is.missing], collapse=", "), sep="")
        if (missing=="fail") {
            stop(msg)
        } else if (missing=="warn") {
            warning(msg)
        }
    }

    return(unlist(unname(nd)))

})
