#'
#' Concatenate genes into a single matrix
#'
#' These functions concatenate separate DNA alignments into a single alignement matrix.
#' \code{concatenate} is a generic with methods for:
#' \itemize{
#'   \item \code{multidna}: returns a \code{DNAbin} matrix
#'   \item \code{multiphyDat}: returns a \code{phyDat} object
#' }
#'
#' @docType methods
#'
#' @rdname concatenate
#'
#' @aliases concatenate
#' @aliases concatenate.generic
#' @aliases concatenate.multidna
#' @aliases concatenate.multiphyDat
#' @aliases concatenate,multidna-method
#' @aliases concatenate,multiphyDat-method
#'
#' @param x a \linkS4class{multidna} or a \linkS4class{multiphyDat} object.
#' @param ... further arguments passed to other methods (currently not used).
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
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
#' image(concatenate(x))
setGeneric("concatenate", function(x, ...) standardGeneric("concatenate"))

#' @rdname concatenate
#'
#' @export
#'
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna is acceptable; by default, all genes are used.
setMethod("concatenate", "multidna", function(x, genes=TRUE, ...){
    x <- add.gaps(x)
    out <- do.call(cbind.DNAbin, x@dna[genes])
    return(out[x@labels,,drop=FALSE])
})

#' @rdname concatenate
#'
#' @export
#'
setMethod("concatenate", "multiphyDat", function(x, genes=TRUE, ...){
    x <- add.gaps(x)
    out <- do.call(cbind.phyDat, x@seq[genes])
    return(out)
})
