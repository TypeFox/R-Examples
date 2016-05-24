
######################
####  CONSTRUCTOR ####
######################
#' multiphyDat constructor
#'
#' New \linkS4class{multiphyDat} objects can be created using \code{new("multiphyDat", ...)} where "..." are arguments documented below.
#' The main input is a list of phyDat matrices. The constructor ensures that all matrices will be reordered in the same way, and genes with missing individuals will be filled by sequences of gaps ("-").
#'
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
#' @aliases initialize,multiphyDat-methods new.multiphyDat
#'
#' @param .Object the object skeleton, automatically generated when calling \code{new}.
#' @param seq a list of phyDat matrices (1 per gene); rows should be labelled and indicate individuals, but different individuals and different orders can be used in different matrices.
#' @param type a character string indicating the type of the sequences stored: "DNA" for DNA sequences, "AA" for amino-acids.
#' @param ind.info an optional data.frame containing information on the individuals, where individuals are in rows.
#' @param gene.info an optional data.frame containing information on the genes, where genes are in rows.
#' @param add.gaps a logical indicating if gap-only sequences should be used where sequences are missing; defaults to TRUE.
#' @param quiet a logical indicating if messages should be shown; defaults to FALSE.
#' @param ... further arguments to be passed to other methods
#'
#' @seealso
#' \itemize{
#' \item the \linkS4class{multiphyDat} class
#' \item \code{\link{read.multiphyDat}}
#' }
#' @examples
#' data(Laurasiatherian)
#' #' ## empty object
#' new("multiphyDat")
#'
#' ## simple conversion with nicely ordered output
#' \dontrun{
#' genes <- list(gene1=subset(Laurasiatherian,, 1:1600, FALSE),
#'     gene2=subset(Laurasiatherian,, 1601:3179, FALSE))
#' x <- new("multiphyDat", genes)
#' x
#' }
#'
#' ## trickier conversion with missing sequences / wrong order
#' genes <- list(gene1=subset(Laurasiatherian, 1:40),
#'     gene2=subset(Laurasiatherian, 8:47))
#' x <- new("multiphyDat", genes)
#' x
#'
setMethod("initialize", "multiphyDat", function(.Object, seq=NULL, type=character(0), ind.info=NULL, gene.info=NULL, add.gaps=TRUE, quiet=FALSE, ...) {

    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object


    ## ESCAPE IF NO DATA ##
    if(is.null(seq)) return(x)


    ## HANDLE SEQ ##
    ## cases where an multiphyDat is provided ##
    if(inherits(seq, "multiphyDat")){
        ind.info <- seq@ind.info
        gene.info <- seq@gene.info
        seq <- seq@seq
        seq@seq <- NULL
        invisible(gc())
    }

    ## cases where no info provided ##
    if(is.null(seq)) return(x)
    if(is.matrix(seq)) seq <- list(seq)

    ## coerce items in SEQ to matrices ##
    # seq <- lapply(seq, as.matrix)
    fun <- function(x)ifelse(is.matrix(x),nrow(x),length(x))
    N.SEQ <- sum(sapply(seq, fun))
    if(N.SEQ==0){
        x@seq <- NULL
        x@ind.info <- x@gene.info <- NULL
        return(x)
    }

    ## convert matrices of characters into phyDat ##
    N.GENES <- length(seq)
    for(i in 1:N.GENES){
        if(is.character(seq[[i]])) seq[[i]] <- phyDat(seq[[i]])
    }

    ## replace with generic names if needed ##
    if(is.null(names(seq))) names(seq) <- paste("gene", 1:N.GENES, sep=".")


    ## get list of all labels ##
    all.labels <- unique(unlist(lapply(seq, names)))
    N.IND <- length(all.labels)


    ## PROCESS META INFO ##
    ## ind.info
    if(!is.null(ind.info)){
        if(nrow(ind.info)>N.IND && !quiet) warning("[multiphyDat constructor] ind.info has more rows than there are individuals")
        ind.info <- as.data.frame(ind.info)
    }

    ## gene.info
    if(!is.null(gene.info)){
        if(nrow(gene.info)>N.GENES && !quiet) warning("[multiphyDat constructor] gene.info has more rows than there are genes")
        gene.info <- as.data.frame(gene.info)
    }


     ## FORM FINAL OUTPUT ##
    x@seq <- seq
    x@type <- as.character(type)
    x@labels <- all.labels
    x@n.ind <- N.IND
    x@n.seq <- as.integer(sum(sapply(x@seq, length)))
    x@ind.info <- ind.info
    x@gene.info <- gene.info

    ## ADD GAPS-ONLY SEQUENCES IF NEEDED ##
    if(add.gaps) x <- add.gaps(x)
    x@n.seq.miss <- .nMissingSequences(x@seq)

    return(x)
}) # end multiphyDat constructor
