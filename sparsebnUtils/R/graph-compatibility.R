#
#  graph-compatibility.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 7/26/15.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#
# Various utility functions for enforcing compatibility with the 'graph' package from BioConductor.
#

#' Conversion to graphNEL object
#'
#' \code{to_graphNEL} converts an object to a \code{\link[graph]{graphNEL}} object. Works on both fitted
#' objects and graphs themselves. In the first case, every underlying 'edges' component is converted to
#' \code{\link[graph]{graphNEL}}. In the second, the conversion applies directly to the object.
#'
#' @param x An object of type \code{\link{sparsebnPath}}, \code{\link{sparsebnFit}}, \code{\link{edgeList}},
#' \code{\link[igraph]{igraph}}, or \code{\link[network]{network}}.
#'
#' @export
to_graphNEL <- function(x){
    ### This function require the 'graph' package to be installed
    if (!requireNamespace("graph", quietly = TRUE)) {
        stop("graph package (from BioConductor) required to coerce data to 'graphNEL' type!", call. = FALSE)
    }

    UseMethod("to_graphNEL", x)
}

#' @export
to_graphNEL.graphNEL <- function(x){
    x
}

#' @export
to_graphNEL.edgeList <- function(x){
    el.graphNEL <- edgeList_to_graphNEL_edgeL(x)
    names.graphNEL <- as.character(1:num.nodes(x))

    graph::graphNEL(nodes = names.graphNEL, edgeL = el.graphNEL, edgemode = 'directed')
}

#' @export
to_graphNEL.igraph <- function(x){
    to_graphNEL(to_edgeList(x))
}

#' @export
to_graphNEL.network <- function(x){
    to_graphNEL(to_edgeList(x))
}

#' @export
to_graphNEL.sparsebnFit <- function(x){
    x$edges <- to_graphNEL(x$edges)

    x
}

#' @export
to_graphNEL.sparsebnPath <- function(x){
    sparsebnPath(lapply(x, to_graphNEL))
}

#
# Helper function to convert an edge list to a graphNEL compatible edge list
#  The main difference is instead of listing parents for each node, graphNEL requires
#  listing children for each parent. There are also different naming and indexing conventions.
#
edgeList_to_graphNEL_edgeL <- function(el){
    #----------- EXAMPLE -----------
    # Default:
    # [[1]]
    # integer(0)
    #
    # [[2]]
    # [1] 1
    #
    # [[3]]
    # [1] 1
    #
    # graphNEL:
    # $`1`
    # $`1`$edges
    # [1] 2 3
    #
    #
    # $`2`
    # $`2`$edges
    # NULL
    #
    #
    # $`3`
    # $`3`$edges
    # NULL
    #
    #-------------------------------

    numnode <- num.nodes(el) # Number of nodes should be same as length of default edge list

    ### Invert the child-parent relationships (aka implicit transpose of adjacency matrix)
    el.graphNEL <- vector(mode = "list", length = numnode)
    for(j in 1:numnode){
        this.column <- el[[j]] # "column" = interpret as column in adj matrix
        for(i in seq_along(this.column)){
            el.graphNEL[[this.column[i]]] <- c(el.graphNEL[[this.column[i]]], j)
        }
    }

    ### Needs an extra component called "edges" (allows for possible specification of weights as well)
    el.graphNEL <- lapply(el.graphNEL, function(x){ list(edges = x)})

    ### List names MUST be character 1,...,numnode
    names(el.graphNEL) <- as.character(1:numnode)

    el.graphNEL
}

#' @export
to_edgeList.graph <- function(x){
    edgeList(graphNEL_to_edgeList_list(x))
}

#
# Helper function to convert a graphNEL object to an edgeList compatible list
#
graphNEL_to_edgeList_list <- function(gr){
    ### Use inEdges to invert the parent-child relationships and coerce the indices to integers
    gr.edgeL <- lapply(graph::inEdges(gr), as.integer)

    ### Eliminate node names
    names(gr.edgeL) <- NULL

    gr.edgeL
}
