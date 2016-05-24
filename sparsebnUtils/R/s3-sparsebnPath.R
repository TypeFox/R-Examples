#
#  s3-sparsebnPath.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#------------------------------------------------------------------------------#
# sparsebnPath S3 Class for R
#------------------------------------------------------------------------------#

#
# sparsebnPath S3 class skeleton
#
# Data
# * <wrapper for a list>
#
# Methods
# * is.sparsebnPath
# * sparsebnPath.list
# * print.sparsebnPath
# * as.list.sparsebnPath
# * num.nodes.sparsebnPath
# * num.edges.sparsebnPath
# * num.samples.sparsebnPath
# * lambda.grid.sparsebnPath
# * get.adjacency.matrix.sparsebnPath
#

#' sparsebnPath class
#'
#' Convenience wrapper class for output of CCDr algorithm: Represents the entire solution path
#' of the CCDr algorithm. Its components are of type \code{\link{sparsebnFit}}. Also inherits
#' from \code{\link{list}}.
#'
#' Each value of lambda in the (discrete) solution path corresponds to a single DAG estimate (see \href{http://jmlr.org/papers/v16/aragam15a.html}{Aragam and Zhou (2015)} for details).
#' Internally, this estimate is represented by a \code{\link{sparsebnFit}} object. The full solution
#' path is then represented as a \code{\link{list}} of \code{\link{sparsebnFit}} objects: This class is essentially a wrapper for this list.
#'
#' @param x Only used internally.
#' @param ... (optional) additional arguments.
#'
#' @section Methods:
#' \code{\link{get.adjacency.matrix}}, \code{\link{lambda.grid}},
#' \code{\link{num.nodes}}, \code{\link{num.edges}}, \code{\link{num.samples}}
#'
#' @docType class
#' @name sparsebnPath
NULL

#' @rdname sparsebnPath
#' @export
is.sparsebnPath <- function(x){
    inherits(x, "sparsebnPath")
} # END IS.sparsebnPath

# sparsebnPath constructor
#' @export
sparsebnPath.list <- function(x){
    if(!check_list_class(x, "sparsebnFit")){
        stop("Some component is not of type sparsebnPath -- sparsebnPath objects must consist of sparsebnFit components only.")
    }

    ### Note that we still allow these objects to inherit from the base list class
    structure(x, class = c("sparsebnPath", "list"))
} # END sparsebnPath.LIST

#' @param verbose If \code{TRUE}, then each estimate in the solution path is printed separately. Do not use for
#'        large graphs or large solution paths. (default = \code{FALSE})
#'
#' @rdname sparsebnPath
#' @method print sparsebnPath
#' @export
print.sparsebnPath <- function(x, verbose = FALSE, ...){
    if(verbose){
        print.default(x) # default generic reverts to list => separate calls to print.sparsebnFit for each component
    } else{
        cat("CCDr solution path\n",
            length(x), " estimates for lambda in [", min(lambda.grid(x)), ",", max(lambda.grid(x)), "]\n",
            "Number of edges per solution: ", paste(num.edges(x), collapse = "-"), "\n",
            num.nodes(x), " nodes\n",
            num.samples(x), " observations\n",
            sep = "")
    }
} # END PRINT.SPARSEBNPATH

#' @export
as.list.sparsebnPath <- function(x, ...){
    class(x) <- "list"
    x
} # END AS.LIST.SPARSEBNPATH

#' @describeIn num.nodes Extracts the number of nodes of \link{sparsebnPath} object.
#' @export
num.nodes.sparsebnPath <- function(x){
    unique(unlist(lapply(x, function(z) z$pp)))
} # END NUM.NODES.SPARSEBNPATH

#' @describeIn num.edges Extracts the number of edges of \link{sparsebnPath} object.
#' @export
num.edges.sparsebnPath <- function(x){
    ### unique(.) not needed since different estimates should have different # of edges
    unlist(lapply(x, function(z) z$nedge))
} # END NUM.EDGES.SPARSEBNPATH

#' @describeIn num.samples Extracts the number of samples of \link{sparsebnPath} object.
#' @export
num.samples.sparsebnPath <- function(x){
    unique(unlist(lapply(x, function(z) z$nn)))
} # END NUM.SAMPLES.SPARSEBNPATH

#' Extract regularization path from solution path
#'
#' @describeIn lambda.grid Returns a vector of lambda values defining the solution path of a \code{\link{sparsebnPath}} object.
#' @export
lambda.grid.sparsebnPath <- function(x){
    lambdas <- unlist(lapply(x, function(z){ z$lambda}))
    names(lambdas) <- NULL

    lambdas
} # END LAMBDA.GRID.sparsebnPath

#' @describeIn get.adjacency.matrix Retrieves all \code{edges} slots in the solution path, converts to an adjacency matrix, and returns as a list
#' @export
get.adjacency.matrix.sparsebnPath <- function(x){
    lapply(x, get.adjacency.matrix)
} # END GET.ADJACENCY.MATRIX.sparsebnPath

#' @export
"[.sparsebnPath" <- function(x, i){
    ### Needed to ensure that naive subsetting of path objects as a list
    ###  still returns a path object
    sparsebnPath(as.list(x)[i])
}

#' @rdname plot.sparsebnFit
#' @method plot sparsebnPath
#' @export
plot.sparsebnPath <- function(x, ...){
    par.default <- par()
    par(mfrow = n2mfrow(length(x)), # Automatically choose a sensible grid to use
        mai=rep(0,4)                   # Need to reset margins (why??? graph packages seem to handle this oddly)
        )

    ### Issues when plotting null DAG, so remove it
    x <- x[-1]

    tryCatch({
        # lapply(x, plot)
        for(fit in x){
            plot(fit, ...)
        }
    }, error = function(c){
        dev.off()
        stop(c)
    })

    par(par.default) # restore user's original settings
}

#' @export
to_edgeList.sparsebnPath <- function(x){
    sparsebnPath(lapply(x, to_edgeList))
}

