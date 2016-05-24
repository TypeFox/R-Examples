## ####################
## ####  ACCESSORS ####
## ####################
#' @name accessors
#' @title multidna Accessors
#' @description Accessors for slots in \linkS4class{multidna} and \linkS4class{multiphyDat} objects.
#'
#' @param x a \linkS4class{multidna} or \linkS4class{multiphyDat} object.
#' @param loci a character, numeric, or logical vector identifying which
#'   loci to return.
#' @param ids a character, numeric, or logical vector identifying which
#'   sequences to return within a locus.
#' @param simplify logical. If \code{FALSE}, always return a list of
#'   DNAbin sequences. If \code{TRUE} and only one locus has been requested,
#'   return a single DNAbin object.
#' @param exclude.gap.only logical. Remove or ignore sequences containing all gaps?
#' @param value a replacement value for the slot.
#' @param ... further arguments passed on to other functions.
#'
#' @details
#' \describe{
#'   \item{getNumInd}{Returns the number of individuals.}
#'   \item{getNumLoci}{Returns the number of loci.}
#'   \item{getLocusNames}{Returns the names of each locus.}
#'   \item{setLocusNames}{Sets the names of each locus.}
#'   \item{getNumSequences}{Returns the number of sequences in each locus.}
#'   \item{getSequenceNames}{Returns the names of individual sequences at each
#'     locus.}
#'   \item{getSequences}{Returns sequences of specified loci and individuals.}
#' }
#'
#setClass("multidna")


## ###############
## ## getNumInd ##
## ###############
#' @rdname accessors
#' @export
setGeneric("getNumInd", function(x, ...) standardGeneric("getNumInd"))
#' @rdname accessors
#' @aliases getNumInd
#' @aliases getNumInd,multidna
#' @export
setMethod("getNumInd", "multidna", function(x, ...) {
  if(is.null(x@dna)) return(0)
  return(x@n.ind)
})
#' @rdname accessors
#' @aliases getNumInd,multiphyDat
#' @export
setMethod("getNumInd", "multiphyDat", function(x, ...) {
    if(is.null(x@seq)) return(0)
    return(x@n.ind)
})


## ################
## ## getNumLoci ##
## ################
#' @rdname accessors
#' @export
setGeneric("getNumLoci", function(x, ...) standardGeneric("getNumLoci"))
#' @rdname accessors
#' @aliases getNumLoci,multidna
#' @export
setMethod("getNumLoci", "multidna", function(x, ...) {
  if(is.null(x@dna)) return(0)
  return(length(x@dna))
})
#' @rdname accessors
#' @aliases getNumLoci,multiphyDat
#' @export
setMethod("getNumLoci", "multiphyDat", function(x, ...) {
    if(is.null(x@seq)) return(0)
    return(length(x@seq))
})


## ###################
## ## getLocusNames ##
## ###################
#' @rdname accessors
#' @export
setGeneric("getLocusNames", function(x, ...) standardGeneric("getLocusNames"))
#' @rdname accessors
#' @aliases getLocusNames,multidna
#' @export
setMethod("getLocusNames", "multidna", function(x, ...) names(x@dna))
#' @rdname accessors
#' @aliases getLocusNames,multiphyDat
#' @export
setMethod("getLocusNames", "multiphyDat", function(x, ...) names(x@seq))


## ###################
## ## setLocusNames ##
## ###################
#' @rdname accessors
#' @export
setGeneric("setLocusNames<-", function(x, value) standardGeneric("setLocusNames<-"))
#' @rdname accessors
#' @aliases setLocusNames<-,multidna
#' @export
setMethod("setLocusNames<-", "multidna", function(x, value) {
  names(x@dna) <- value
  validObject(x)
  x
})
#' @rdname accessors
#' @aliases setLocusNames<-,multiphyDat
#' @export
setMethod("setLocusNames<-", "multiphyDat", function(x, value) {
    names(x@seq) <- value
    validObject(x)
    x
})


## #####################
## ## getNumSequences ##
## #####################
#' @rdname accessors
#' @export
setGeneric("getNumSequences", function(x, ...) standardGeneric("getNumSequences"))
#' @rdname accessors
#' @aliases getNumSequences,multidna
#' @export
setMethod("getNumSequences", "multidna",
          function(x, exclude.gap.only = TRUE, loci = NULL, ...) {
  # check that object isn't empty
  if(is.null(x@dna)) {
    warning("'x' is empty. NULL returned.", call. = FALSE)
    return(NULL)
  }

  if(is.character(loci) | is.null(loci)) {
    loci <- .checkLocusNames(x, loci)
  } else if(is.numeric(loci) | is.logical(loci)) {
    loci <- getLocusNames(x)[loci]
  } else stop("'loci' must be a character, numeric, or logical")

  sapply(loci, function(this.locus) {
    dna <- x@dna[[this.locus]]
    if(exclude.gap.only) sum(!.isGapOnly(dna)) else nrow(as.matrix(dna))
  })
})
#' @rdname accessors
#' @aliases getNumSequences,multiphyDat
#' @export
setMethod("getNumSequences", "multiphyDat",
          function(x, exclude.gap.only = TRUE, loci = NULL, ...) {
              # check that object isn't empty
              if(is.null(x@seq)) {
                  warning("'x' is empty. NULL returned.", call. = FALSE)
                  return(NULL)
              }

              if(is.character(loci) | is.null(loci)) {
                  loci <- .checkLocusNames(x, loci)
              } else if(is.numeric(loci) | is.logical(loci)) {
                  loci <- getLocusNames(x)[loci]
              } else stop("'loci' must be a character, numeric, or logical")

              sapply(loci, function(this.locus) {
                  dna <- x@seq[[this.locus]]
                  if(exclude.gap.only) sum(!.isGapOnlyPhyDat(dna)) else length(dna)
              })
          })


## ######################
## ## getSequenceNames ##
## ######################
#' @rdname accessors
#' @export
setGeneric("getSequenceNames", function(x, ...) standardGeneric("getSequenceNames"))
#' @rdname accessors
#' @aliases getSequenceNames,multidna
#' @export
setMethod("getSequenceNames", "multidna",
          function(x, exclude.gap.only = TRUE, loci = NULL, ...) {
  # check that object isn't empty
  if(is.null(x@dna)) {
    warning("'x' is empty. NULL returned.", call. = FALSE)
    return(NULL)
  }

  if(is.character(loci) | is.null(loci)) {
    loci <- .checkLocusNames(x, loci)
  } else if(is.numeric(loci) | is.logical(loci)) {
    loci <- getLocusNames(x)[loci]
  } else stop("'loci' must be a character, numeric, or logical")

  sapply(loci, function(this.locus) {
    dna <- x@dna[[this.locus]]
    if(exclude.gap.only) labels(dna)[!.isGapOnly(dna)] else labels(dna)
  }, simplify = FALSE)
})
#' @rdname accessors
#' @aliases getSequenceNames,multiphyDat
#' @export
setMethod("getSequenceNames", "multiphyDat",
          function(x, exclude.gap.only = TRUE, loci = NULL, ...) {
              # check that object isn't empty
              if(is.null(x@seq)) {
                  warning("'x' is empty. NULL returned.", call. = FALSE)
                  return(NULL)
              }

              if(is.character(loci) | is.null(loci)) {
                  loci <- .checkLocusNames(x, loci)
              } else if(is.numeric(loci) | is.logical(loci)) {
                  loci <- getLocusNames(x)[loci]
              } else stop("'loci' must be a character, numeric, or logical")

              sapply(loci, function(this.locus) {
                  dna <- x@seq[[this.locus]]
                  if(exclude.gap.only) labels(dna)[!.isGapOnlyPhyDat(dna)] else labels(dna)
              }, simplify = FALSE)
          })


## ##################
## ## getSequences ##
## ##################
#' @rdname accessors
#' @export
setGeneric("getSequences", function(x, ...) standardGeneric("getSequences"))
#' @rdname accessors
#' @aliases getSequences,multidna
#' @export
setMethod("getSequences", "multidna",
          function(x, loci = NULL, ids = NULL, simplify = TRUE,
                   exclude.gap.only = TRUE, ...) {
  # check that object isn't empty
  if(is.null(x@dna)) {
    warning("'x' is empty. NULL returned.", call. = FALSE)
    return(NULL)
  }

  if(is.character(loci) | is.null(loci)) {
    loci <- .checkLocusNames(x, loci)
  } else if(is.numeric(loci) | is.logical(loci)) {
    loci <- getLocusNames(x)[loci]
  } else stop("'loci' must be a character, numeric, or logical")

  # loop through loci
  new.dna <- sapply(loci, function(this.locus) {
    # extract this DNAbin object
    dna <- as.list(x@dna[[this.locus]])
    if(exclude.gap.only) dna <- dna[!.isGapOnly(dna)]
    # return sequences for IDs which are present
    locus.ids <- .checkIDs(dna, ids)
    if(is.null(locus.ids)) NULL else dna[locus.ids]
  }, simplify = FALSE)
  new.dna <- new.dna[!sapply(new.dna, is.null)]

  if(length(new.dna) == 1 & simplify) new.dna[[1]] else new.dna
})
#' @rdname accessors
#' @aliases getSequences,multiphyDat
#' @export
setMethod("getSequences", "multiphyDat",
          function(x, loci = NULL, ids = NULL, simplify = TRUE,
                   exclude.gap.only = TRUE, ...) {
              # check that object isn't empty
              if(is.null(x@seq)) {
                  warning("'x' is empty. NULL returned.", call. = FALSE)
                  return(NULL)
              }

              if(is.character(loci) | is.null(loci)) {
                  loci <- .checkLocusNames(x, loci)
              } else if(is.numeric(loci) | is.logical(loci)) {
                  loci <- getLocusNames(x)[loci]
              } else stop("'loci' must be a character, numeric, or logical")
              # loop through loci
              new.dna <- sapply(loci, function(this.locus) {
                  # extract this DNAbin object
                  dna <- as.list(x@seq[[this.locus]])
                  if(exclude.gap.only) dna <- subset(dna, !.isGapOnlyPhyDat(dna))
                  # return sequences for IDs which are present
                  locus.ids <- .checkIDs(dna, ids)
                  if(is.null(locus.ids)) NULL else subset(dna,locus.ids)
              }, simplify = FALSE)
              new.dna <- new.dna[!sapply(new.dna, is.null)]
              if(length(new.dna) == 1 & simplify) new.dna[[1]] else new.dna
          })