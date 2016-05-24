#'
#' Remove gap-only sequences for missing data
#'
#' In \linkS4class{multidna} and \linkS4class{multiphyDat}, some individuals may not be sequenced for all genes, resulting in gap-only sequences for missing data.
#' The generic function \code{rm.gaps} has method for both objects; it identifies the missing sequences, and removes gap-only sequences from the alignments wherever needed.
#'
#'
#' @docType methods
#'
#' @export
#'
#' @aliases rm.gaps
#'
#' @param x a \linkS4class{multidna} or \linkS4class{multiphyDat} object.
#' @param ... further arguments passed to other methods (currently not used).
#'
#' @aliases rm.gaps.generic
#' @aliases rm.gaps.multidna
#' @aliases rm.gaps.multiphyDat
#' @aliases rm.gaps,multidna-method
#' @aliases rm.gaps,multiphyDat-method
#'
setGeneric("rm.gaps", function(x, ...) standardGeneric("rm.gaps"))



#' @rdname rm.gaps
#'
#' @export
#'
setMethod("rm.gaps", "multidna", function(x, ...){
    ## ESCAPE IF NO DNA SEQUENCES ##
    if(is.null(x@dna)) return(x)

    ## FUNCTION TO REORDER MATRIX AND REMOVE MISSING SEQUENCES ##
    form.dna.matrix <- function(mat.dna, labels){
        ## find sequences to remove
        toRemove <- apply(mat.dna==as.DNAbin("-"),1,all)

        ## return relevant sequences
        return(mat.dna[!toRemove, , drop=FALSE])
    }

    ## APPLY THIS FUNCTION TO ALL MATRICES ##
    x@dna <- lapply(x@dna, form.dna.matrix, x@labels)

    ## REMOVE EMPTY MATRICES ##
    x@dna <- x@dna[sapply(x@dna, nrow)>0]

    ## UPDATE NUMBER OF SEQUENCES ##
    x@n.seq <- as.integer(sum(sapply(x@dna, nrow)))
    x@n.seq.miss <- .nMissingSequences(x@dna)

    ## RETURN OBJECT ##
    return(x)
}) # end multidna method




#' @rdname rm.gaps
#'
#' @export
#'
setMethod("rm.gaps", "multiphyDat", function(x, ...){
    ## ESCAPE IF NO DNA SEQUENCES ##
    if(is.null(x@seq)) return(x)

    ## FUNCTION TO REORDER MATRIX AND REMOVE MISSING SEQUENCES ##
    form.dna.phyDat <- function(dna, labels){
        ## convert dna to matrix of characters
        dna <- as.character(dna)

        ## find sequences to remove
        toRemove <- apply(dna=="-",1,all)

        ## return relevant sequences
        return(as.phyDat(dna[!toRemove,,drop=FALSE]))
    }

    ## APPLY THIS FUNCTION TO ALL MATRICES ##
    x@seq <- lapply(x@seq, form.dna.phyDat, x@labels)

    ## REMOVE EMPTY OBJECTS ##
    x@seq <- x@seq[sapply(x@seq, length)>0]

    ## UPDATE NUMBER OF SEQUENCES ##
    x@n.seq <- as.integer(sum(sapply(x@seq, length)))
    x@n.seq.miss <- .nMissingSequences(x@seq)

    ## RETURN OBJECT ##
    return(x)
})
