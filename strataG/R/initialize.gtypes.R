#' @title \code{gtypes} Constructor
#' @description Create a new \linkS4class{gtypes} object using 
#'   \code{new("gtypes", ...)}, where '\code{...}' are arguments 
#'   documented below.
#' 
#' @param .Object the object skeleton, automatically generated when 
#'   calling \code{new}.
#' @param gen.data a vector, matrix, or data.frame containing the alleles 
#'   at each locus. See below for more details.
#' @param ploidy ploidy of the loci.
#' @param ind.names an optional vector of individual sample names.
#' @param sequences an optional \linkS4class{multidna} 
#'   object containing sequences represented by each locus.
#' @param strata an optional stratification scheme from \code{schemes}.
#' @param schemes an optional data.frame of stratification schemes.
#' @param description an optional description for the object.
#' @param other other optional information to include.
#' @param remove.sequences logical. If \code{TRUE} any sequences not referenced 
#'
#' @details
#' For multi-allele loci, the \code{gen.data} argument should be 
#' formatted such that every consecutive \code{ploidy} columns represent 
#' alleles of one locus. Locus names are taken from the column names in 
#' \code{gen.data} and should be formatted with the same root locus name, with 
#' unique suffixes representing alleles (e.g., for Locus1234: Locus1234.1 
#' and Locus1234.2, or Locus1234_A and Locus1234_B). \cr\cr
#' If \code{gen.data} is a vector it is assumed to represent haplotypes of a 
#' haploid marker.
#' Sample names can be either in the rownames of \code{gen.data} or given 
#' separately in \code{ind.names}. If \code{ind.names} are provided, these are 
#' used in lieu of rownames in \code{gen.data} and also used to label 
#' \code{schemes}.
#' If sequences are provided in \code{sequences}, then they should be named 
#' and match haplotype labels in \code{loc.col} of \code{x}. If multiple 
#' genes are given as a \linkS4class{multidna}, then they should have the 
#' same names as column names in \code{X} from \code{loc.col} to the end.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{df2gtypes}, \link{sequence2gtypes}, \link{gtypes2df},
#'   \link{gtypes2genind}, \link{gtypes2loci}
#'
#' @aliases initialize.gtypes new
#' @importFrom methods setMethod
#' 
setMethod("initialize", "gtypes", 
          function(.Object, gen.data, ploidy, ind.names = NULL,
                   sequences = NULL, strata = NULL, schemes = NULL, 
                   description = NULL, other = NULL, 
                   remove.sequences = FALSE) {
  
  if(is.null(gen.data) | is.null(ploidy)) return(.Object)
            
  # check gen.data
  if(is.vector(gen.data)) {
    gen.data <- cbind(gen.data)
    colnames(gen.data) <- "Haplotype"
  }
  if(!(is.matrix(gen.data) | is.data.frame(gen.data))) {
    stop("'gen.data' is not a vector, matrix, or data.frame")
  }
  
  # check ploidy
  ploidy <- as.integer(ploidy)
  if(ncol(gen.data) %% ploidy != 0) {
    stop("the number of columns in 'gen.data' is not a multiple of 'ploidy'")
  }
  if(ploidy > 1 & !is.null(sequences)) {
    stop("'sequences' can't be present if 'ploidy' is greater than 1")
  }
  
  # check ind.names
  ind.names <- if(!is.null(ind.names)) {
    if(!is.vector(ind.names)) stop("'ind.names' must be a vector")
    if(length(ind.names) != nrow(gen.data)) {
      stop("the length of 'ind.names' must equal the number of rows in 'gen.data'")
    }
    as.character(ind.names)
  } else {
    if(!is.null(rownames(gen.data))) rownames(gen.data) else 1:nrow(gen.data)
  }
  if(any(duplicated(ind.names))) {
    dup.names <- unique(ind.names[duplicated(ind.names)])
    dup.names <- paste(dup.names, collapse = ", ")
    stop(paste("there are duplicated individual names:", dup.names))
  }
  rownames(gen.data) <- ind.names
  
  # check schemes
  if(!is.null(schemes)) {
    if(!(is.matrix(schemes) | is.data.frame(schemes))) {
      stop("'schemes' is not a matrix or data.frame")
    } else {
      schemes <- as.data.frame(schemes)
      for(i in 1:ncol(schemes)) schemes[, i] <- factor(schemes[, i])
      if(is.null(rownames(schemes))) {
        if(nrow(schemes) == length(ind.names)) {
          rownames(schemes) <- ind.names
        } else {
          stop("'schemes' doesn't have rownames and not as long as 'ind.names'")
        }
      } else {
        if(length(intersect(rownames(schemes), rownames(gen.data))) == 0) {
          stop("no rownames in 'schemes' are in 'gen.data'")
        }
      }
    }
  }
  rm(ind.names)
  
  # check description
  if(is.null(description)) {
    description <- paste(
      "gtypes created on", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  }
  
  # format loci
  locus.cols <- matrix(1:ncol(gen.data), nrow = ploidy)
  loci <- do.call(data.frame, lapply(1:ncol(locus.cols), function(i) {
    factor(unlist(gen.data[, locus.cols[, i]], use.names = FALSE))
  }))
  colnames(loci) <- if(is.null(colnames(gen.data))) {
    # return generic names if no colnames assigned
    nums <- formatC(1:ncol(loci), digits = floor(log10(ncol(loci))), flag = "0")
    paste("Locus", nums, sep = "_")
  } else .parseLocusNames(colnames(gen.data), ploidy)
  row.ids <- rep(rownames(gen.data), ploidy)
  row.nums <- rep(1:ploidy, each = nrow(gen.data))
  rownames(loci) <- paste(row.ids, row.nums, sep = ".")
  rm(gen.data, locus.cols, row.ids, row.nums)
  
  # create valid locus names
  colnames(loci) <- make.names(colnames(loci))
  
  # check sequences
  if(!is.null(sequences)) {
    sequences <- as.multidna(sequences)
    if(getNumLoci(sequences) != ncol(loci)) {
      stop("the number of genes in 'sequences' is not equal to the number of loci")
    }
    setLocusNames(sequences) <- colnames(loci)
    for(loc in colnames(loci)) {
      haps <- unique(as.character(loci[, loc]))
      seq.names <- getSequenceNames(sequences)[[loc]]
      missing <- setdiff(haps, seq.names)
      if(length(missing) > 0) {
        missing <- paste(missing, collapse = ", ")
        stop(paste("the following haplotypes can't be found in sequences for locus '", loc, "': ", missing, sep = ""))
      }
    }
  }
  
  # create and return gtypes object
  g <- .Object
  g@loci <- loci
  g@ploidy <- ploidy
  g@sequences <- sequences
  g@schemes <- schemes
  g@description <- description
  g@other <- other
  
  # Check for samples missing data for all loci
  g.mat <- as.matrix(g)
  not.missing.all <- apply(g.mat, 1, function(y) !all(is.na(y)))
  to.keep <- names(which(not.missing.all))
  to.remove <- names(which(!not.missing.all))
  if(length(to.remove)) {
    warning("The following samples are missing data for all loci and have been removed: ", 
            paste(to.remove, collapse = ", "))
  }
  g@loci <- g@loci[idRows(to.keep, rownames(g@loci)), , drop = FALSE]
  g@loci <- droplevels(g@loci)
  
  if(remove.sequences) g <- removeSequences(g)
  
  stratify(g, strata)
})
         