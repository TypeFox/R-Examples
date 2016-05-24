#' Pairwise distances for multiple gene data
#'
#' This function computes pairwise genetic distances between individuals using genes in a \linkS4class{multidna} object.
#' By default, one distance matrix (dist object) is created for each each, but a single distance can be derived by pooling all genes (argument \code{pool=TRUE})
#'
#' @export
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a \linkS4class{multidna} object.
#' @param pool a logical indicating if all genes should be pooled (concatenated) to obtain a single distance matrix; defaults to FALSE.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna is acceptable; by default, all genes are used.
#' @param ... further arguments passed to \code{\link[ape]{dist.dna}}.
#'
#' @aliases dist.multidna
#'
#' @seealso \code{\link[ape]{dist.dna}}
#'
#' @return a list of dist objects (pool=FALSE) or a single dist object (pool=TRUE)
#'
#' @examples
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' plot(x)
#'
#' ## get separate distance matrix and pooled one
#' lD <- dist.multidna(x)
#' D <- dist.multidna(x, pool=TRUE)
#'
#' ## get corresponding NJ trees
#' ltrees <- lapply(lD, nj)
#' tree <- nj(D)
#'
#' par(mfrow=c(3,1))
#' for(i in 1:2) plot(ltrees[[i]], main=names(ltrees)[i])
#' plot(tree, main="Pooled distances")
#'
dist.multidna <- function(x, pool=FALSE, genes=TRUE, ...){
    ## subset data
    x <- x[,genes]

    ## if genes are  pooled
    if(pool){
        x <- concatenate(x)
        out <- dist.dna(x, ...)
        return(out)
    }

    ## otherwise: one tree per gene
    out <- lapply(x@dna, dist.dna, ...)
    names(out) <- names(x@dna)
    return(out)
}
