
######################
####  CONSTRUCTOR ####
######################
#' multidna constructor
#'
#' New \linkS4class{multidna} objects can be created using \code{new("multidna", ...)} where "..." are arguments documented below.
#' The main input is a list of DNAbin matrices. The constructor ensures that all matrices will be reordered in the same way, and as an option (setting \code{add.gaps=TRUE}, gap-only sequences ("...-----...") will be added wherever sequences are missing.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
#' @aliases initialize,multidna-methods new.multidna
#'
#' @param .Object the object skeleton, automatically generated when calling \code{new}.
#' @param dna a list of DNAbin matrices (1 per gene); rows should be labelled and indicate individuals, but different individuals and different orders can be used in different matrices.
#' @param ind.info an optional data.frame containing information on the individuals, where individuals are in rows.
#' @param gene.info an optional data.frame containing information on the genes, where genes are in rows.
#' @param add.gaps a logical indicating if gap-only sequences should be used where sequences are missing; defaults to TRUE.
#' @param quiet a logical indicating if messages should be shown; defaults to FALSE.
#' @param ... further arguments to be passed to other methods
#'
#' @seealso
#' \itemize{
#' \item the \linkS4class{multidna} class
#' \item \code{\link{read.multidna}} and \code{\link{read.multidna}}
#' }
#' @examples
#'
#' ## empty object
#' new("multidna")
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' image(woodmouse)
#' image(x@@dna[[1]])
#' image(x@@dna[[2]])
#'
#' ## trickier conversion with missing sequences / wrong order
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
#' x <- new("multidna", genes)
#' x
#' image(x@@dna[[1]])
#' image(x@@dna[[2]])
#'
setMethod("initialize", "multidna", function(.Object, dna=NULL, ind.info=NULL, gene.info=NULL, add.gaps=TRUE, quiet=FALSE, ...) {

    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object


    ## ESCAPE IF NO DATA ##
    if(is.null(dna)) return(x)


    ## HANDLE DNA ##
    ## cases where an multidna is provided ##
    if(inherits(dna, "multidna")){
        ind.info <- dna@ind.info
        gene.info <- dna@gene.info
        dna <- dna@dna
        dna@dna <- NULL
        invisible(gc())
    }

    ## cases where no info provided ##
    if(is.null(dna)) return(x)
    if(is.matrix(dna)) dna <- list(dna)

    ## coerce items in DNA to matrices ##
    dna <- lapply(dna, as.matrix)
    N.SEQ <- sum(sapply(dna, nrow))
    if(N.SEQ==0){
        x@dna <- NULL
        x@ind.info <- x@gene.info <- NULL
        return(x)
    }

    ## convert matrices of characters into DNAbin ##
    N.GENES <- length(dna)
    for(i in 1:N.GENES){
        if(is.character(dna[[i]])) dna[[i]] <- as.DNAbin(dna[[i]])
    }

    ## replace with generic names if needed ##
    if(is.null(names(dna))) names(dna) <- paste("gene", 1:N.GENES, sep=".")


    ## HANDLE LABELS ##
    ## handle missing labels ##
    missing.labels <- any(sapply(dna, function(e) is.null(labels(e))))
    if(missing.labels){
        if(!quiet) message("[multidna constructor] missing/incomplete labels provided - using generic labels.\n")
        ## error if varying numbers of rows
        if(length(unique(sapply(dna, nrow)))>1) stop("[multidna constructor] no labels provided and varying number of sequences across genes - cannot assume individuals are identical.")
        labels <- paste("individual", 1:nrow(dna[[1]]), sep=".")
        for(i in 1:N.GENES) rownames(dna[[i]]) <- labels
    }

    ## get list of all labels ##
    all.labels <- unique(unlist(lapply(dna,rownames)))
    N.IND <- length(all.labels)


    ## SORT MATRICES OF DNA ##
    sort.mat <- function(mat){
        if(is.null(rownames(mat))) return(mat)
        return(mat[sort(rownames(mat)),,drop=FALSE])
    }
    dna <- lapply(dna, sort.mat)


    ## PROCESS META INFO ##
    ## ind.info
    if(!is.null(ind.info)){
        if(nrow(ind.info)>N.IND && !quiet) warning("[multidna constructor] ind.info has more rows than there are individuals")
        ind.info <- as.data.frame(ind.info)
    }

    ## gene.info
    if(!is.null(gene.info)){
        if(nrow(gene.info)>N.GENES && !quiet) warning("[multidna constructor] gene.info has more rows than there are genes")
        gene.info <- as.data.frame(gene.info)
    }


    ## FORM FINAL OUTPUT ##
    x@dna <- dna
    x@labels <- all.labels
    x@n.ind <- N.IND
    x@n.seq <- as.integer(sum(sapply(x@dna, nrow)))
    x@ind.info <- ind.info
    x@gene.info <- gene.info

    ## ADD GAPS-ONLY SEQUENCES IF NEEDED ##
    if(add.gaps) x <- add.gaps(x)
    x@n.seq.miss <- .nMissingSequences(x@dna)

    return(x)
}) # end multidna constructor
