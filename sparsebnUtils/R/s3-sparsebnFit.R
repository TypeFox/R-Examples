#
#  s3-sparsebnFit.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#------------------------------------------------------------------------------#
# sparsebnFit S3 Class for R
#------------------------------------------------------------------------------#

#
# sparsebnFit S3 class skeleton
#
# Data
# * edgeList edges          // edge list, adjacency matrix, or graphNEL object of DAG estimate
# * numeric lambda          // regularization parameter
# * integer nedge           // number of edges
# * integer pp              // number of nodes
# * integer nn              // number of observations
# * numeric time            // time to run CCDr algorithm
#
# Methods
# * is.sparsebnFit
# * sparsebnFit.list
# * as.list.sparsebnFit
# * print.sparsebnFit
# * get.adjacency.matrix
# * num.nodes.sparsebnFit
# * num.edges.sparsebnFit
# * num.samples.sparsebnFit
# * to_B.sparsebnFit
#

#' sparsebnFit class
#'
#' Main class for representing DAG estimates: Represents a single DAG estimate in the solution path.
#' Generally speaking, these estimates should be wrapped up in a \code{\link{sparsebnPath}} object, but
#' can be handled separately if desired (be careful!).
#'
#' @param x Only used internally.
#' @param ... (optional) additional arguments.
#'
#' @section Slots:
#' \describe{
#' \item{\code{edges}}{(edgeList) Edge list of estimated DAG (see \code{\link{edgeList}}).}
#' \item{\code{lambda}}{(numeric) Value of lambda for this estimate.}
#' \item{\code{nedge}}{(integer) Number of edges in this estimate.}
#' \item{\code{pp}}{(integer) Number of nodes.}
#' \item{\code{nn}}{(integer) Number of observations this estimate was based on.}
#' \item{\code{time}}{(numeric) Time in seconds to generate this estimate.}
#' }
#'
#' @section Methods:
#' \code{\link{get.adjacency.matrix}}
#' \code{\link{num.nodes}}, \code{\link{num.edges}}, \code{\link{num.samples}}
#'
#' @docType class
#' @name sparsebnFit
NULL

#' @rdname sparsebnFit
#' @export
is.sparsebnFit <- function(x){
    inherits(x, "sparsebnFit")
} # END IS.sparsebnFit

# sparsebnFit constructor
#' @method sparsebnFit list
#' @export
sparsebnFit.list <- function(x){

    #
    # Need to be careful when using this constructor directly since it allows the nedge
    #  component to be different from the actual number of edges stored in the SBM object.
    #  This is allowed for efficiency reasons while running the main algorithm.
    #
    # UPDATE: An explicit check has been added for now.
    #

    if( !is.list(x)){
        stop("Input must be a list!")
    } else if( length(x) != 6 || !setequal(names(x), c("edges", "lambda", "nedge", "pp", "nn", "time"))){
        stop("Input is not coercable to an object of type sparsebnFit, check list for the following elements: edges (edgeList), lambda (numeric), nedge (integer), pp (integer), nn (integer), time (numeric or NA)")
    } else if( !is.edgeList(x$edges)){
        stop("'edges' component must be a valid edgeList object!")
    } else if(num.edges(x$edges) != x$nedge){
        stop("Attempting to set nedge to an improper value: Must be equal to the number of nonzero values in edges.")
    }

    ### Update values to be consistent with edgeList
    if(x$pp != num.nodes(x$edges)){
        stop("Attempting to create sparsebnFit object with inconsistent number of nodes! input = ", x$pp, " != output = ", num.nodes(x$edges))
    }
    x$pp <- num.nodes(x$edges)

    if(x$nedge != num.edges(x$edges)){
        stop("Attempting to create sparsebnFit object with inconsistent number of edges! input = ", x$nedge, " != output = ", num.edges(x$edges))
    }
    x$nedge <- num.edges(x$edges)

    ### Output with DAG as edgeList
    out <- structure(x, class = "sparsebnFit")

    ### Coerce to user's desired data structure
    pkg_graph <- getGraphPackage()
    if(!is.null(pkg_graph)){
        if(pkg_graph == "graph"){
            out <- to_graphNEL(out)
        } else if(pkg_graph == "igraph"){
            out <- to_igraph(out)
        } else if(pkg_graph == "network"){
            out <- to_network(out)
        } else{
            stop(invalid_pkg_specification())
        }
    }

    out
} # END sparsebnFit.LIST

#' @method as.list sparsebnFit
#' @export
as.list.sparsebnFit <- function(x, ...){
    list(edges = x$edges, lambda = x$lambda, nedge = x$nedge, pp = x$pp, nn = x$nn, time = x$time)
} # END AS.LIST.sparsebnFit
#' @rdname sparsebnFit
#' @method print sparsebnFit
#' @export
print.sparsebnFit <- function(x, ...){
    MAX_NODES <- 20

    cat("CCDr estimate\n",
        x$nn, " observations\n",
        "lambda = ", x$lambda, "\n",
        sep = "")

    cat("\nDAG: \n")
    print(x$edges)
    if(x$pp < MAX_NODES) {
        # print(get.adjacency.matrix(x))
    }
} # END PRINT.sparsebnFit

#' @describeIn get.adjacency.matrix Retrieves \code{edges} slot and converts to an adjacency matrix
#' @export
get.adjacency.matrix.sparsebnFit <- function(x){
    get.adjacency.matrix.edgeList(x$edges)
} # END GET.ADJACENCY.MATRIX.sparsebnFit

#' @describeIn num.nodes Extracts the number of nodes of \link{sparsebnFit} object.
#' @export
num.nodes.sparsebnFit <- function(x){
    x$pp
} # END NUM.NODES.sparsebnFit

#' @describeIn num.edges Extracts the number of edges of \link{sparsebnFit} object.
#' @export
num.edges.sparsebnFit <- function(x){
    x$nedge
} # END NUM.EDGES.sparsebnFit

#' @describeIn num.samples Extracts the number of samples of \link{sparsebnFit} object.
#' @export
num.samples.sparsebnFit <- function(x){
    x$nn
} # END NUM.SAMPLES.sparsebnFit

#' Plot a fitted Bayesian network object
#'
#' Plot plot plot
#'
#' @param x fitted object (\code{\link{sparsebnFit}} or \code{\link{sparsebnPath}}) to plot.
#' @param ... (optional) additional arguments to plotting mechanism.
#'
#' @method plot sparsebnFit
#' @export
plot.sparsebnFit <- function(x, ...){
    pkg_plot <- getPlotPackage()

    if(!is.null(pkg_plot)){
        if(pkg_plot == "graph"){
            plot(to_graphNEL(x$edges), ...)
        } else if(pkg_plot == "igraph"){
            plot(to_igraph(x$edges), ...)
        } else if(pkg_plot == "network"){
            plot(to_network(x$edges), ...)
        } else{
            stop("Incorrect package specified. Must be one of: 'graph', 'igraph', 'network'.")
        }
    } else{
        stop("No package specified for plotting! This is an internal error and should not happen -- please report this issue.")
    }
}

#' @describeIn to_edgeList description
#' @export
to_edgeList.sparsebnFit <- function(x){
    x$edges <- to_edgeList(x$edges)

    x
}

#------------------------------------------------------------------------------#
# to_B.sparsebnFit
# Internal function to convert estimates from the (Rho, R) parametrization to
#  the standard (B, Omega) parametrization.
#
# !!! 1-29-16: This function needs to be deprecated
#
# to_B.sparsebnFit <- function(fit){
#     .Deprecated()
#     fit$sbm <- to_B(fit$sbm)
#
#     fit
# }
