#
#  network-compatibility.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/30/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#
# Various utility functions for enforcing compatibility with the 'network' package (part of statnet).
#

#' Conversion to network object
#'
#' \code{to_network} converts an object to a \code{\link[network]{network}} object. Works on both fitted
#' objects and graphs themselves. In the first case, every underlying 'edges' component is converted to
#' \code{\link[network]{network}}. In the second, the conversion applies directly to the object.
#'
#' @param x An object of type \code{\link{sparsebnPath}}, \code{\link{sparsebnFit}}, \code{\link{edgeList}},
#' \code{\link[graph]{graphNEL-class}}, or \code{\link[igraph]{igraph}}
#'
#' @export
to_network <- function(x){
    ### This function requires the 'network' package to be installed
    if (!requireNamespace("network", quietly = TRUE)) {
        stop("network package required to coerce data to 'network' type!", call. = FALSE)
    }

    UseMethod("to_network", x)
}

#' @export
to_network.network <- function(x){
    x
}

#' @export
to_network.edgeList <- function(x){
    if(num.edges(x) > 0){
        el.network <- edgeList_to_network_edgelist(x)

        network::network(el.network, directed = TRUE)
    } else{
        #
        # Is there a way to intialize a null graph using the network constructor?
        #
        network::network.initialize(num.nodes(x)) # hack to work around special case of null graph
    }

}

#' @export
to_network.graphNEL <- function(x){
    to_network(to_edgeList(x))
}

#' @export
to_network.igraph <- function(x){
    to_network(to_edgeList(x))
}

#' @export
to_network.sparsebnFit <- function(x){
    x$edges <- to_network(x$edges)

    x
}

#' @export
to_network.sparsebnPath <- function(x){
    sparsebnPath(lapply(x, to_network))
}

#
# Helper function to convert an edge list to a network compatible edge list matrix
#  The output is a 2 column matrix where each row represents one edge:
#       [i j] = (i, j)
#
# NOTE: The required input is the _same_ as for igraph, so we just hijack that function.
#
edgeList_to_network_edgelist <- function(el){
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
    # network:
    # 1 2
    # 1 3
    #
    #-------------------------------

    edgeList_to_igraph_edgelist(el)
}

#' @export
to_edgeList.network <- function(x){
    edgeList(network_to_edgeList_list(x))
}

#
# Helper function to convert a network object to an edgeList compatible list
#
network_to_edgeList_list <- function(net){
    net.edgelist <- network::as.edgelist(net)
    numnode <- network::network.size(net)

    edgelist_mat_to_edgeList_list(net.edgelist, numnode)
### Deprecated in favour of edgelist_mat_to_edgeList_list
#     net.edgeL <- lapply(vector("list", length = numnode), as.integer)
#     if(network::network.edgecount(net) > 0){ # If no edges, simply return list full of integer(0)'s
#         for(j in 1:nrow(net.edgelist)){
#             ### NOTE: Fix this to be memory-efficient for large graphs
#             net.edgeL[[net.edgelist[j, 2]]] <- c(net.edgeL[[net.edgelist[j, 2]]], net.edgelist[j, 1])
#         }
#     }
#
#     net.edgeL
}
