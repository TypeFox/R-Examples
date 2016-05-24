#
#  s3-edgeList.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#------------------------------------------------------------------------------#
# edgeList S3 Class for R
#------------------------------------------------------------------------------#

#
# edgeList S3 class skeleton
#
# Data
# * <wrapper for a list>
#
# Methods
# * is.edgeList
# * edgeList.list
# * print.edgeList
# * as.matrix.edgeList
# * as.list.edgeList
# * get.adjacency.matrix.edgeList
# * num.nodes.edgeList
# * num.edges.edgeList
# * is.zero.edgeList
#

#' edgeList class
#'
#' Convenience wrapper class for a (column-major) edge list. Each component of the list
#' corresponds to a node, and each component is an integer vector whose components are the parents
#' of this node in the graph.
#'
#' Also inherits from \code{\link{list}}.
#'
#' @param x A list containing parents for each node in a graph. The length of this list
#'          should be the same as the number of nodes in the graph.
#'
#' @section Methods:
#' \code{\link{get.adjacency.matrix}},
#' \code{\link{num.nodes}}, \code{\link{num.edges}}
#'
#' @docType class
#' @name edgeList
NULL

#' @rdname edgeList
#' @export
is.edgeList <- function(x){
    inherits(x, "edgeList")
} # END IS.EDGELIST

#' @export
edgeList.list <- function(x){
    if(!is.list(x)){
        stop("Input must be a list!")
    }

    ### Don't allow length zero lists: What would that represent anyway?
    if(length(x) == 0){
        stop("Input must have at least one component!")
    }

    ### Cannot assign a parent larger than total number of nodes, or < 1
    max.node.index <- suppressWarnings(max(unlist(x))) # Ignore warning if graph is empty
    min.node.index <- suppressWarnings(min(unlist(x))) #
    if(max.node.index > length(x) || min.node.index < 1){
        stop(sprintf("The indices of the parents must be between 1 and number of nodes (= length of list)!
                      Currently between %d and %d.", min.node.index, max.node.index))
    }

    structure(x, class = c("edgeList", "list"))
} # END EDGELIST.LIST

#' @method print edgeList
#' @export
print.edgeList <- function(x, ...){
    if(num.edges(x) == 0){
        edgeL.out <- paste0("<Empty graph on ", num.nodes(x), " nodes.>")
    } else{
        ### Assumes the DAG has at most 1000 nodes: Output will be cramped and illegible if the graph is larger than this
        ### We shouldn't be printing this when pp > 1000 anyway!
        edgeL.out <- mapply(function(x, y){
            prefix <- paste0("[", x, "]")
            prefix <- sprintf("%-5s", prefix)
            paste0(prefix, paste(sprintf("%4d", sort(y)), collapse = ""))
        }, 1L:length(x), x)
        edgeL.out <- unlist(edgeL.out)
        edgeL.out <- paste(edgeL.out, collapse = " \n")
    }

    cat("edgeList object\n", edgeL.out, "\n", sep = "")
}

#' @export
as.matrix.edgeList <- function(x, ...){
    as.matrix(get.adjacency.matrix.edgeList(x))
}

#' @export
as.list.edgeList <- function(x, ...){
    class(x) <- "list"
    x
}

#' @describeIn get.adjacency.matrix Convert internal \code{edgeList} representation to an adjacency matrix
#' @export
get.adjacency.matrix.edgeList <- function(x){
    numnode <- length(x)
    Matrix.out <- Matrix::Matrix(0, nrow = numnode, ncol = numnode)

    ### This loop is pretty slow!
    for(j in seq_along(x)){
        for(i in x[[j]]){
            Matrix.out[i, j] <- 1
        }
    }

    Matrix.out
} # END GET.ADJACENCY.MATRIX.EDGELIST

#' @describeIn num.nodes Extracts the number of nodes of \link{edgeList} object.
#' @export
num.nodes.edgeList <- function(x){
    length(x)
} # END NUM.NODES.EDGELIST

#' @describeIn num.edges Extracts the number of edges of \link{edgeList} object.
#' @export
num.edges.edgeList <- function(x){
    sum(sapply(x, length))
} # END NUM.EDGES.EDGELIST

#' @describeIn is.zero Determines whether or not the object represents a null graph with no edges.
#' @export
is.zero.edgeList <- function(x){
    (num.edges(x) == 0)
} # END IS.ZERO.EDGELIST

#' @export
to_edgeList.edgeList <- function(x){
    x
}

#
# Convert a standard two-column edge list to an edgeList compatible list
#
edgelist_mat_to_edgeList_list <- function(x, numnode){

    edgeL <- lapply(vector("list", length = numnode), as.integer)
    if(nrow(x) > 0){ # If no edges, simply return list full of integer(0)'s
        for(j in 1:nrow(x)){
            ### NOTE: Fix this to be memory-efficient for large graphs
            edgeL[[x[j, 2]]] <- c(edgeL[[x[j, 2]]], x[j, 1])
        }
    }

    edgeL
}

