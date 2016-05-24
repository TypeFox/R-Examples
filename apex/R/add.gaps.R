#'
#' Add gap-only sequences for missing data
#'
#' In \linkS4class{multidna} and \linkS4class{multiphyDat}, some individuals may not be sequenced for all genes.
#' The generic function \code{add.gaps} has method for both objects; it identifies the missing sequences, and adds gap-only sequences to the alignments wherever needed.
#'
#'
#' @docType methods
#'
#' @export
#'
#' @aliases add.gaps
#'
#' @param x a \linkS4class{multidna} or \linkS4class{multiphyDat} object.
#' @param ... further arguments passed to other methods (currently not used).
#'
#' @aliases add.gaps.generic
#' @aliases add.gaps.multidna
#' @aliases add.gaps.multiphyDat
#' @aliases add.gaps,multidna-method
#' @aliases add.gaps,multiphyDat-method
#'
setGeneric("add.gaps", function(x, ...) standardGeneric("add.gaps"))



#' @rdname add.gaps
#'
#' @export
#'
setMethod("add.gaps", "multidna", function(x, ...){
    ## ESCAPE IF NO DNA SEQUENCES ##
    if(is.null(x@dna)) return(x)

    ## FUNCTION TO REORDER MATRIX AND ADD MISSING SEQUENCES ##
    form.dna.matrix <- function(mat.dna, labels){
        ## make matrix of missing sequences if needed
        lab.missing <- labels[!(labels %in% rownames(mat.dna))]
        n.missing <- length(lab.missing)
        if(n.missing>0){
            mat.NA <- as.DNAbin(matrix("-", ncol=ncol(mat.dna), nrow=n.missing))
            rownames(mat.NA) <- lab.missing
            mat.dna <- rbind(mat.dna, mat.NA)
        }

        ## return ordered sequences
        return(mat.dna[labels,])
    }

    ## APPLY THIS FUNCTION TO ALL MATRICES ##
    x@dna <- lapply(x@dna, form.dna.matrix, x@labels)

    ## update number of sequences ##
    x@n.seq <- as.integer(sum(sapply(x@dna, nrow)))
    x@n.seq.miss <- .nMissingSequences(x@dna)

    ## RETURN OBJECT ##
    return(x)
}) # end multidna method




#' @rdname add.gaps
#'
#' @export
#'
setMethod("add.gaps", "multiphyDat", function(x, ...){
    ## ESCAPE IF NO DNA SEQUENCES ##
    if(is.null(x@seq)) return(x)

    ## FUNCTION TO REORDER MATRIX AND ADD MISSING SEQUENCES ##
    form.dna.phyDat <- function(dna, labels){
        ## convert dna to matrix of characters
        dna <- as.character(dna)

        ## make matrix of missing sequences if needed
        lab.missing <- labels[!(labels %in% labels(dna))]
        n.missing <- length(lab.missing)
        if(n.missing>0){
            mat.NA <- matrix("-", ncol=ncol(dna), nrow=n.missing)
            rownames(mat.NA) <- lab.missing
            dna <- rbind(dna, mat.NA)
        }

        ## return ordered sequences
        return(as.phyDat(dna[labels,]))
    }

    ## APPLY THIS FUNCTION TO ALL MATRICES ##
    x@seq <- lapply(x@seq, form.dna.phyDat, x@labels)

    ## update number of sequences ##
    x@n.seq <- as.integer(sum(sapply(x@seq, length)))
    x@n.seq.miss <- .nMissingSequences(x@seq)

    ## RETURN OBJECT ##
    return(x)
})
