
######################
####  SHOW METHOD ####
######################
#' Display multidna objects
#'
#' Default printing for multidna objects
#'
#' @export
#'
#' @aliases show,multidna-method
#' @aliases show.multidna
#'
#' @param object a multidna object
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
setMethod ("show", "multidna", function(object){
    ## INFO FOR PRINTING ##
    x <- object
    n.genes <- ifelse(is.null(x@dna), 0, length(x@dna))
    seqword <- ifelse(x@n.seq>1, "sequences", "sequence")
    seqmissword <- ifelse(x@n.seq.miss>1, "sequences", "sequence")
    geneword <- ifelse(n.genes>1, "genes", "gene")
    indword <-  ifelse(x@n.ind>1, "individuals", "individual")

    ## PRINT OBJECT ##
    cat(paste("=== multidna ===\n",sep=""))
    cat(paste("[",x@n.seq,"DNA", seqword, "in", n.genes, geneword,"]\n"))

    cat("\n@n.ind:", x@n.ind, indword)
    cat("\n@n.seq:", x@n.seq, seqword, "in total")
    cat("\n@n.seq.miss:", x@n.seq.miss, "gap-only (missing)", seqmissword)
    cat("\n@labels:", head(x@labels))
    if(length(x@labels>6)) cat("...")

    cat("\n")
    if(n.genes>0) {
        cat("\n@dna: (list of DNAbin matrices)\n")
        print(object@dna)
    }

    if(!is.null(x@ind.info)) {
        cat("\n@ind.info:", nrow(x@ind.info), "rows,", ncol(x@ind.info), "columns\n")
        print(head(x@ind.info))
        if(nrow(x@ind.info)>6) cat("...\n")
    }
    if(!is.null(x@gene.info)) {
        cat("\n@gene.info:", nrow(x@gene.info), "rows,", ncol(x@gene.info), "columns\n")
        print(head(x@gene.info))
        if(nrow(x@gene.info)>6) cat("...")
    }

    cat("\n")
})







######################
####  SHOW METHOD ####
######################
#' Display multiphyDat objects
#'
#' Default printing for multiphyDat objects
#'
#' @export
#'
#' @aliases show,multiphyDat-method
#' @aliases show.multiphyDat
#'
#' @param object a multiphyDat object
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
setMethod ("show", "multiphyDat", function(object){
    ## INFO FOR PRINTING ##
    x <- object
    n.genes <- ifelse(is.null(x@seq), 0, length(x@seq))
    seqword <- ifelse(x@n.seq>1, "sequences", "sequence")
    seqmissword <- ifelse(x@n.seq.miss>1, "sequences", "sequence")
    geneword <- ifelse(n.genes>1, "genes", "gene")
    indword <-  ifelse(x@n.ind>1, "individuals", "individual")

    ## PRINT OBJECT ##
    cat(paste("=== multiphyDat ===\n",sep=""))
    cat(paste("[",x@n.seq,"DNA", seqword, "in", n.genes, geneword,"]\n"))
    cat("\n@type:", x@type)
    cat("\n@n.ind:", x@n.ind, indword)
    cat("\n@n.seq:", x@n.seq, seqword, "in total")
    cat("\n@n.seq.miss:", x@n.seq.miss, "gap-only (missing)", seqmissword)
    cat("\n@labels:", head(x@labels))
    if(length(x@labels>6)) cat("...")

    cat("\n")
    if(n.genes>0) {
        cat("\n@seq: (list of phyDat objects)\n")
        print(object@seq)
    }

    if(!is.null(x@ind.info)) {
        cat("\n@ind.info:", nrow(x@ind.info), "rows,", ncol(x@ind.info), "columns\n")
        print(head(x@ind.info))
        if(nrow(x@ind.info)>6) cat("...\n")
    }
    if(!is.null(x@gene.info)) {
        cat("\n@gene.info:", nrow(x@gene.info), "rows,", ncol(x@gene.info), "columns\n")
        print(head(x@gene.info))
        if(nrow(x@gene.info)>6) cat("...")
    }

    cat("\n")
})





