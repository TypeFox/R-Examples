
## TODO -- the behavior of edgeLength needs to be made more consistent
## with other functions like MRCA. The user should be able to specify a
## vector of nodes, of edges, or both.

##### This file contains
## hasEdgeLength
## edgeLength and edgeLength<-
## isUltrametric
## nodeDepth
## sumEdgeLength

##' edgeLength methods
##'
##' These functions give information about and allow replacement of edge lengths.
##'
##' The \code{edgeLength} function returns the edge length in the same
##' order as the edges in the matrix.
##'
##' @param x a \code{phylo4} or \code{phylo4d} object.
##' @param value a numeric vector indicating the new values for the edge lengths
##' @param node optional numeric or character vector indicating the
##' nodes for which edge
##' @param use.names should the the name attributes of \code{value} be
##' used to match the length to a given edge.
##' @param tol the tolerance to decide whether all the tips have the
##' same depth to test if the tree is ultrametric. Default is
##' \code{.Machine$double.eps^0.5}.
##' @param \dots optional arguments (none used at present).
##' @return \describe{
##'
##' \item{hasEdgeLength}{whether or not the object has edge lengths
##' (logical)}
##'
##' \item{edgeLength}{a named vector of the edge length for the
##' object}
##'
##' \item{isUltrametric}{whether or not the tree is ultrametric (all
##' the tips are have the same depth (distance from the root) (logical)}
##'
##' \item{sumEdgeLength}{the sum of the edge lengths for a set of
##' nodes (intended to be used with \code{ancestors} or \code{descendants})}
##'
##' \item{nodeHeight}{the distance between a node and the root or the
##' tips. The format of the result will depend on the options and the
##' number of nodes provided, either a vector or a list.}
##'
##'  \item{nodeDepth}{Deprecated, now replaced by \code{nodeHeight}. A
##' named vector indicating the \dQuote{depth} (the distance between
##' the root and a given node).}
##'
##' \item{depthTip}{Deprecated, now  replaced by \code{nodeHeight}.}
##'
##' }
##' @seealso \code{ancestors}, \code{descendants}, \code{.Machine} for
##' more information about tolerance.
##' @export
##' @docType methods
##' @aliases hasEdgeLength
##' @rdname edgeLength-methods
##' @include phylo4-class.R
##' @include phylo4-methods.R
##' @include nodeId-methods.R
##' @examples
##'   data(geospiza)
##'   hasEdgeLength(geospiza) # TRUE
##'   topoGeo <- geospiza
##'   edgeLength(topoGeo) <- NULL
##'   hasEdgeLength(topoGeo)  # FALSE
##'
##'   edgeLength(geospiza)[2]       # use the position in vector
##'   edgeLength(geospiza)["16-17"] # or the name of the edge
##'   edgeLength(geospiza, 17)      # or the descendant node of the edge
##'
##'   ## The same methods can be used to update an edge length
##'   edgeLength(geospiza)[2] <- 0.33
##'   edgeLength(geospiza)["16-17"] <- 0.34
##'   edgeLength(geospiza, 17) <- 0.35
##'
##'   ## Test if tree is ultrametric
##'   isUltrametric(geospiza)   # TRUE
##'   ## indeed all tips are at the same distance from the root
##'   nodeHeight(geospiza, nodeId(geospiza, "tip"), from="root")
##'   ## compare distances from tips of two MRCA
##'   nodeHeight(geospiza, MRCA(geospiza, c("pallida", "psittacula")), from="min_tip")
##'   nodeHeight(geospiza, MRCA(geospiza, c("fortis", "difficilis")), from="min_tip")
##'   ## or the same but from the root
##'   nodeHeight(geospiza, MRCA(geospiza, c("pallida", "psittacula")), from="root")
##'   nodeHeight(geospiza, MRCA(geospiza, c("fortis", "difficilis")), from="root")
setGeneric("hasEdgeLength", function(x) {
    standardGeneric("hasEdgeLength")
})

##' @rdname edgeLength-methods
##' @aliases hasEdgeLength,phylo4-method
setMethod("hasEdgeLength", signature(x="phylo4"),
 function(x) {
    !all(is.na(x@edge.length))
})

#### ----- edgeLength

##' @rdname edgeLength-methods
##' @aliases edgeLength
##' @export
setGeneric("edgeLength", function(x, ...) {
    standardGeneric("edgeLength")
})

##' @rdname edgeLength-methods
##' @aliases edgeLength,phylo4-method
setMethod("edgeLength", signature(x="phylo4"),
 function(x, node) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    if (!missing(node)) {
        id <- getEdge(x, node)
    } else {
        id <- edgeId(x, "all")
    }
    elen <- x@edge.length[match(id, names(x@edge.length))]
    names(elen) <- id
    return(elen)
})

##' @rdname edgeLength-methods
##' @aliases edgeLength<-
##' @export
setGeneric("edgeLength<-", function(x, use.names=TRUE, ..., value) {
    standardGeneric("edgeLength<-")
})

##' @name edgeLength<-
##' @rdname edgeLength-methods
##' @aliases edgeLength<-,phylo4-method edgeLength<-,phylo4,ANY-method
setReplaceMethod("edgeLength", signature(x="phylo4", value="ANY"),
 function(x, use.names=TRUE, ..., value) {
    len <- .createEdge(value, x@edge, type="lengths", use.names)
    ## return empty vector if all values are NA
    if (all(is.na(len))) {
        emptyVec <- numeric()
        attributes(emptyVec) <- list(names=character(0))
        x@edge.length <- emptyVec
    } else {
        x@edge.length <- len
    }
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
})

##### ------ depthTips

##' @rdname edgeLength-methods
##' @aliases depthTips
##' @export
setGeneric("depthTips", function(x) {
  standardGeneric("depthTips")
})

##' @rdname edgeLength-methods
##' @aliases depthTips,phylo4-methods
setMethod("depthTips", signature(x="phylo4"), function(x) {
              .Deprecated("nodeHeight")
  nodeDepth(x, 1:nTips(x))
})

##### ----- nodeDepth

##' @rdname edgeLength-methods
##' @aliases nodeDepth
##' @export
setGeneric("nodeDepth", function(x, node) {
  standardGeneric("nodeDepth")
})

##' @rdname edgeLength-methods
##' @aliases nodeDepth,phylo4-method
setMethod("nodeDepth", signature(x="phylo4"),
          function(x, node) {
              .Deprecated("nodeHeight")
    if (!hasEdgeLength(x))
      return(NULL)
    else {
      node <- getNode(x, node, missing="warn")
      node <- node[!is.na(node)]
      res <- sapply(node, function(n)
                    sumEdgeLength(x, ancestors(x, n, "ALL")))
      if (length(res) == 1) {
        res <- res[[1]]
        names(res) <- names(node)
      }
      res
    }
})

###### ----- nodeHeight

##' @param from The point of reference for calculating the height of
##' the node. \code{root} calculates the distance between the root of
##' the tree and the node. \code{all_tip} return the distance between
##' the node and all the tips descending from it. \code{min_tip} the
##' distance between the node and its closest tip. \code{max_tip} the
##' distance between the node and its farther tip. \code{min_tip} and
##' \code{max_tip} will be identical if the tree is ultrametric. If
##' more than one tip is equidistant from the node, the tip with the
##' lowest node id will be returned.
##' @rdname edgeLength-methods
##' @aliases nodeHeight
##' @export
setGeneric("nodeHeight", function(x, node, from) {
               standardGeneric("nodeHeight")
           })

##' @rdname edgeLength-methods
##' @aliases nodeHeight,phylo4-method
setMethod("nodeHeight", signature(x = "phylo4"),
          function(x, node, from = c("root", "all_tip", "min_tip", "max_tip")) {
              from <- match.arg(from)
              if (!hasEdgeLength(x))
                  return(NULL)
              else {
                  node <- getNode(x, node, missing = "warn")
                  node <- node[!is.na(node)]
                  tip_id <- nodeId(x, "tip")
                  if (from != "root") {
                      ## Get the full paths to the tips from the node
                      res <- lapply(node, function(n) {
                                        if (n %in% tip_id) {
                                            ## tips are always at 0
                                            tmp_res <- stats::setNames(0, tipLabels(x)[n])
                                        } else {
                                            desc_pths <- descendants(x, n, "all")
                                            ## Get the paths in the other direction
                                            anc_pths <- lapply(desc_pths[desc_pths %in% tip_id],
                                                               function(n) {
                                                                   ancestors(x, n, "ALL")
                                                               })
                                            ## Shortest paths for each tip
                                            pths <- lapply(anc_pths, function(anc_pth) {
                                                               intersect(desc_pths, anc_pth)
                                                           })
                                            tmp_res <- sapply(pths, function(n) {
                                                                  sumEdgeLength(x, n)
                                                              })

                                            tmp_res <- switch(from,
                                                              "all_tip" = tmp_res,
                                                              "min_tip" = tmp_res[which.min(tmp_res)],
                                                              "max_tip" = tmp_res[which.max(tmp_res)])
                                        }
                                        tmp_res
                                    })
                  } else {
                      res <- sapply(node, function(n) {
                                        sumEdgeLength(x, ancestors(x, n, "ALL"))
                                    })
                  }
                  names(res) <- node
                  if (length(res) == 1) {
                      res <- res[[1]]
                  }
              }
              res
          })


###### ----- sumEdgeLength

##' @rdname edgeLength-methods
##' @aliases sumEdgeLength
##' @export
setGeneric("sumEdgeLength", function(x, node) {
    standardGeneric("sumEdgeLength")
})

##' @rdname edgeLength-methods
##' @aliases sumEdgeLength,phylo4-method
setMethod("sumEdgeLength", signature(x="phylo4"),
 function(x, node) {
    if(!hasEdgeLength(x))
        NULL
    else {
        nd <- getNode(x, node)
        iEdges <- which(x@edge[,2] %in% nd)
        sumEdges <- sum(x@edge.length[iEdges], na.rm=TRUE)
        sumEdges
    }
})

###### ----- isUltrametric

##' @rdname edgeLength-methods
##' @aliases isUltrametric
##' @export
setGeneric("isUltrametric", function(x, tol=.Machine$double.eps^.5) {
  standardGeneric("isUltrametric")
})

##' @rdname edgeLength-methods
##' @aliases isUltrametric,phylo4-method
setMethod("isUltrametric", signature(x="phylo4"),
  function(x, tol=.Machine$double.eps^.5) {
    if (!hasEdgeLength(x)) {
      stop("The tree has no edge lengths.")
    }
    if (identical(all.equal.numeric(stats::var(nodeHeight(x, nodeId(x, "tip"), "root")), 0, tolerance=tol), TRUE)) {
      TRUE
    }
    else FALSE
  })
