#' @name gtypes.accessors
#' @title \code{gtypes} Accessors
#' @description Accessors for slots in \linkS4class{gtypes} objects.
#' 
#' @param x a \linkS4class{gtypes} object.
#' @param seqName the name (or number) of a set of sequences from the 
#'   \code{@@sequences} slot to return.
#' @param ids vector of individual ids.
#' @param loci vector of loci.
#' @param as.haplotypes return sequences as haplotypes? If \code{TRUE}, contents of 
#'   \code{@@sequences} slot are returned. If \code{FALSE}, one sequence per individual 
#'   is returned.
#' @param i,j,k subsetting slots for individuals (\code{i}), loci (\code{j}),
#'   or strata (\code{k}). See Details for more information.
#' @param quiet suppress warnings about unmatched requested individuals, loci, or strata?
#' @param drop if \code{TRUE} the return object will have unused sequences removed.
#' @param ... other arguments passed from generics (ignored).
#' @param value value being assigned by accessor.
#' 
#' @details 
#' Indexing a \code{gtypes} object with integers, characters, or logicals with 
#'   the \code{[} operator follows the same rules as normal R indexing. The 
#'   order that individuals, loci, and strata are chosen in follow the order 
#'   returned by \code{indNames}, \code{locNames}, and \code{strataNames} 
#'   respectively. If unstratified samples are present, they can be selected as
#'   a group either by including \code{NA} in the character or numeric vector of the 
#'   \code{k} slot, or by providing a logical vector based on \code{is.na(strata(g))} 
#'   to the \code{i} slot.
#'
#' @return
#' \tabular{ll}{
#'   \code{nInd} \tab number of individuals.\cr
#'   \code{nLoc} \tab number of loci.\cr
#'   \code{nStrata} \tab number of strata.\cr
#'   \code{indNames} \tab vector of individual/sample names.\cr
#'   \code{locNames} \tab vector of locus names.\cr
#'   \code{strataNames} \tab vector of strata names for current scheme.\cr
#'   \code{ploidy} \tab number of alleles at each locus.\cr
#'   \code{other} \tab contents of \code{@@other} slot.\cr
#'   \code{strata} \tab return or modify the current stratification.\cr
#'   \code{schemes} \tab return or modify the current stratification schemes.\cr
#'   \code{loci} \tab return a data.frame of the alleles for the specified ids 
#'     and loci.\cr
#'   \code{sequences} \tab return the \linkS4class{multidna} object in 
#'     the \code{@@sequences} slot.\cr See \code{\link[apex]{getSequences}} to 
#'     extract individual genes or sequences from this object.
#'   \code{description} \tab return the object's description.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @aliases accessors
#' 
#' @examples
#' #--- create a diploid (microsatellite) gtypes object
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' nStrata(msats.g)
#' strataNames(msats.g)
#' nLoc(msats.g)
#' locNames(msats.g)
#' 
#' # reassign all samples to two randomly chosen strata
#' strata(msats.g) <- sample(c("A", "B"), nInd(msats.g), rep = TRUE)
#' msats.g
#' 
#' 
#' #--- a sequence example
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' wood.g <- sequence2gtypes(x)
#' strata(wood.g) <- sample(c("A", "B"), nInd(wood.g), rep = TRUE)
#' wood.g
#' 
#' # get the multidna sequence object
#' multi.seqs <- sequences(wood.g)
#' class(multi.seqs) # "multidna"
#'
#' # get a list of DNAbin objects
#' dnabin.list <- getSequences(multi.seqs)
#' class(dnabin.list) # "list"
#' 
#' # get a DNAbin object of the first locus
#' dnabin.1 <- getSequences(multi.seqs, locNames(wood.g)[1])
#' class(dnabin.1) # "DNAbin"
#' 
#' # NOTE: The default to the 'simplify' argument in 'getSequences' is TRUE, 
#' #   so if there is only one locus, 'getSequences' will return a DNAbin object
#' #   rather than a single element list unless 'simplify = FALSE':
#' gene1 <- wood.g[, "gene1", ]
#' gene1.dnabin <- getSequences(sequences(gene1))
#' class(gene1.dnabin) # "DNAbin"
#' 
setClass("gtypes")


#' @rdname gtypes.accessors
#' @aliases nInd
#' @export
#' 
setMethod("nInd", "gtypes", function(x, ...) nrow(x@loci) / x@ploidy)


#' @rdname gtypes.accessors
#' @aliases nLoc
#' @export
#' 
setMethod("nLoc", "gtypes", function(x, ...) ncol(x@loci))


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("nStrata", function(x, ...) standardGeneric("nStrata"))

#' @rdname gtypes.accessors
#' @aliases nStrata
#' @export
#' 
setMethod("nStrata", "gtypes", function(x, ...) nlevels(x@strata))


#' @rdname gtypes.accessors
#' @aliases indNames
#' @export
#' 
setMethod("indNames", "gtypes", function(x, ...) {
  ids <- rownames(x@loci)[1:(nrow(x@loci) / x@ploidy)]
  unique(substr(ids, 1, nchar(ids) - 2))
})


#' @rdname gtypes.accessors
#' @aliases locNames
#' @export
#' 
setMethod("locNames", "gtypes", function(x, ...) colnames(x@loci))


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("strataNames", function(x, ...) standardGeneric("strataNames"))

#' @rdname gtypes.accessors
#' @aliases strataNames
#' @export
#' 
setMethod("strataNames", "gtypes", function(x, ...) levels(x@strata))


#' @rdname gtypes.accessors
#' @aliases ploidy
#' @export
#' 
setMethod("ploidy", "gtypes", function(x, ...) x@ploidy)


#' @rdname gtypes.accessors
#' @aliases other
#' @export
#' 
setMethod("other", "gtypes", function(x, ...) x@other)


#' @rdname gtypes.accessors
#' @aliases strata
#' @export
#' 
setMethod("strata", "gtypes", function(x) x@strata)


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("strata<-", function(x, value) standardGeneric("strata<-"))

#' @rdname gtypes.accessors
#' @aliases strata
#' @importFrom methods validObject
#' @export
#' 
setMethod("strata<-", "gtypes", function(x, value) {
  strata <- factor(rep(value, length.out = nInd(x)))
  names(strata) <- indNames(x)
  x@strata <- droplevels(strata)
  validObject(x)
  x
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("schemes", function(x, ...) standardGeneric("schemes"))

#' @rdname gtypes.accessors
#' @aliases schemes
#' @export
#' 
setMethod("schemes", "gtypes", function(x, ...) x@schemes)


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("schemes<-", function(x, value) standardGeneric("schemes<-"))

#' @rdname gtypes.accessors
#' @aliases schemes
#' @importFrom methods validObject
#' @export
#' 
setMethod("schemes<-", "gtypes", function(x, value) {
  x@schemes <- value
  validObject(x)
  x
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("loci", function(x, ...) standardGeneric("loci"))

#' @rdname gtypes.accessors
#' @aliases loci
#' @export
#' 
setMethod("loci", "gtypes", function(x, ids = NULL, loci = NULL) {
  if(is.null(ids)) ids <- indNames(x)
  if(is.null(loci)) loci <- locNames(x)
  if(!all(ids %in% indNames(x))) stop("some 'ids' not found in 'x'")
  if(!all(loci %in% locNames(x))) stop("some 'loci' not found in 'x'")
  x@loci[idRows(ids, rownames(x@loci)), loci, drop = FALSE]
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("sequences", function(x, ...) standardGeneric("sequences"))

#' @rdname gtypes.accessors
#' @aliases sequences
#' @export
#' 
setMethod("sequences", "gtypes", function(x, seqName = NULL, as.haplotypes = TRUE, ...) {
  if(is.null(x@sequences)) return(NULL)
  dna <- getSequences(x@sequences, simplify = FALSE)
  if(!as.haplotypes) {
    dna <- lapply(locNames(x), function(l) {
      haps <- as.character(loci(x)[, l])
      ind.seqs <- dna[[l]][haps]
      names(ind.seqs) <- indNames(x)
      ind.seqs
    })
  }
  if(!is.null(seqName)) dna <- dna[seqName]
  as.multidna(dna)
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("description", function(x, ...) standardGeneric("description"))

#' @rdname gtypes.accessors
#' @aliases description
#' @export
#' 
setMethod("description", "gtypes", function(x, ...) x@description)


#' @rdname gtypes.accessors
#' @aliases index subset
#' @importFrom methods new
#' @export
#' 
setMethod("[", 
          signature(x = "gtypes", i = "ANY", j = "ANY", drop = "ANY"), 
          function(x, i, j, k, ..., quiet = TRUE, drop = FALSE) {
  
  # check ids (i)
  if(missing(i)) i <- TRUE
  if(is.factor(i)) i <- as.character(i)
  ids <- indNames(x)
  i <- if(is.character(i)) {
    i <- unique(i)
    missing.ids <- setdiff(i, ids)
    if(length(missing.ids) > 0) {
      missing.ids <- paste(missing.ids, collapse = ", ")
      warning(paste("the following ids cannot be found:", missing.ids))
    }
    intersect(i, ids)
  } else ids[i]
  
  # check loci (j)
  if(missing(j)) j <- TRUE
  if(is.factor(j)) j <- as.character(j)
  locs <- locNames(x)
  j <- if(is.character(j)) {
    j <- unique(j)
    missing.locs <- setdiff(j, locs)
    if(length(missing.locs) > 0 & !quiet) {
      missing.locs <- paste(missing.locs, collapse = ", ")
      warning(paste("the following loci cannot be found:", missing.locs))
    }
    intersect(j, locs)
  } else locs[j]
  
  # check strata (k) 
  if(missing(k)) k <- TRUE
  if(is.factor(k)) k <- as.character(k)
  st <- strata(x)
  k.i <- if(is.character(k)) {
    k <- unique(k)
    missing.strata <- setdiff(k, st)
    if(length(missing.strata) > 0 & !quiet) {
      missing.strata <- paste(missing.strata, collapse = ", ")
      warning(paste("the following strata cannot be found:", missing.strata))
    }
    st.i <- names(st)[st %in% k]
    if(any(is.na(k))) st.i <- c(st.i, names(st)[is.na(st)])
    st.i
  } else {
    names(st)[st %in% strataNames(x)[k]]
  }
  
  # check ids in selected strata
  missing.ids <- setdiff(i, k.i)
  if(length(missing.ids) > 0 & !quiet) {
    missing.ids <- paste(missing.ids, collapse = ", ")
    warning(paste("the following ids are not in the selected strata:", missing.ids))
  }
  i <- intersect(i, k.i)
  
  if(length(i) == 0) stop("no samples selected")
  if(length(j) == 0) stop("no loci selected")
  if(length(k) == 0) stop("no strata selected")
  
  i <- i[order(match(i, ids))]
  x@loci <- x@loci[idRows(i, rownames(x@loci)), j, drop = FALSE]
  x@loci <- droplevels(x@loci)
  x@strata <- droplevels(x@strata[i])
  if(!is.null(x@sequences)) {
    j.seqs <- getSequences(x@sequences, loci = j, 
                           simplify = FALSE, exclude.gap.only = FALSE)
    x@sequences <- new("multidna", j.seqs)
  }
  
  # Check for samples missing data for all loci
  x.mat <- as.matrix(x)
  not.missing.all <- apply(x.mat, 1, function(y) !all(is.na(y)))
  to.keep <- names(which(not.missing.all))
  to.remove <- names(which(!not.missing.all))
  if(length(to.remove)) {
    warning("The following samples are missing data for all loci and have been removed: ", 
            paste(to.remove, collapse = ", "))
  }
  x@loci <- x@loci[idRows(to.keep, rownames(x@loci)), , drop = FALSE]
  x@loci <- droplevels(x@loci)
  x@strata <- droplevels(x@strata[to.keep])
  
  if(drop) x <- removeSequences(x)

  return(x)
})
