#
#  igraph-compatibility.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/30/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#
# Various utility functions for enforcing compatibility with the 'igraph' package.
#

#' Conversion to igraph object
#'
#' \code{to_igraph} converts an object to a \code{\link[igraph]{igraph}} object. Works on both fitted
#' objects and graphs themselves. In the first case, every underlying 'edges' component is converted to
#' \code{\link[igraph]{igraph}}. In the second, the conversion applies directly to the object.
#'
#' @param x An object of type \code{\link{sparsebnPath}}, \code{\link{sparsebnFit}}, \code{\link{edgeList}},
#' \code{\link[graph]{graphNEL-class}}, or \code{\link[network]{network}}.
#'
#' @export
to_igraph <- function(x){
    ### This function requires the 'igraph' package to be installed
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("igraph package required to coerce data to 'igraph' type!", call. = FALSE)
    }

    UseMethod("to_igraph", x)
}

#' @export
to_igraph.igraph <- function(x){
    x
}

#' @export
to_igraph.edgeList <- function(x){
    if(num.edges(x) > 0){
        el.igraph <- edgeList_to_igraph_edgelist(x)
        igraph::graph_from_edgelist(el.igraph, directed = TRUE)
    } else{
        ### Special case to create empty graph
        igraph::graph.empty(n = num.nodes(x), directed = TRUE)
    }
}

#' @export
to_igraph.graphNEL <- function(x){
    to_igraph(to_edgeList(x))
}

#' @export
to_igraph.network <- function(x){
    to_igraph(to_edgeList(x))
}

#' @export
to_igraph.sparsebnFit <- function(x){
    x$edges <- to_igraph(x$edges)

    x
}

#' @export
to_igraph.sparsebnPath <- function(x){
    sparsebnPath(lapply(x, to_igraph))
}

#
# Helper function to convert an edge list to a igraph compatible edge list matrix
#  The output is a 2 column matrix where each row represents one edge:
#       [i j] = (i, j)
#
edgeList_to_igraph_edgelist <- function(el){
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
    # igraph:
    # 1 2
    # 1 3
    #
    #-------------------------------

    ### Slick version
    el.igraph <- matrix(NA, nrow=0, ncol=2) # base case for null graph
    if(num.edges(el) > 0){
        #
        # Bugfix: By default, SIMPLIFY is set to TRUE, which causes problem when
        #         do.call is called below ('Error in do.call("rbind", el.igraph) : second argument must be a list')
        #         Need to set SIMPLIFY = FALSE to prevent pre-mature matrixification.
        #
        el.igraph <- mapply(function(x, y){
            if(length(x) > 0) cbind(x, y)
        }, el, 1:num.nodes(el), SIMPLIFY = FALSE)

        el.igraph <- do.call("rbind", el.igraph)
    }

    ### Equivalent code: What the slick version does under the hood
#         numnode <- num.nodes(el) # number of nodes in edgeList
#         numedge <- num.edges(el) # number of edges in edgeList
#
#         el.igraph <- matrix(NA, ncol = 2, nrow = numedge)
#         start <- 1
#         for(j in 1:numnode){
#             this.size <- length(el[[j]])
#             end <- start + this.size - 1
#             if(this.size > 0){
#                 el.igraph[start:end, ] <- cbind(el[[j]], j) # el[[j]] are the parents
#                 start <- end + 1
#             }
#         }

    el.igraph
}

#' @export
to_edgeList.igraph <- function(x){
    edgeList(igraph_to_edgeList_list(x))
}

#
# Helper function to convert a igraph object to an edgeList compatible list
#
igraph_to_edgeList_list <- function(igr){
    igr.edgeL <- igraph::as_adj_list(igr, mode = "in")
    igr.edgeL <- lapply(igr.edgeL, as.integer)

    igr.edgeL
}
