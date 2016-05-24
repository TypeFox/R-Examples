setClassUnion("dnaSequences", c("multidna", "NULL"))

#' @title \code{gtypes} Class
#' @description An S4 class storing multi-allelic locus or sequence data along
#'   with a current stratification and option stratification schemes.
#'
#' @slot loci a data.frame containing the allelic data as one column per locus.
#'   Alleles are on multiple rows per column with samples listed in the same
#'   order for each allele. rownames are sample names plus allele number
#'   formatted as 12345.1 and 12345.2 where 12345 is the sample name and 1 and
#'   2 are the first and second alleles. colnames are unique locus names.
#' @slot sequences a \linkS4class{multidna} object.
#' @slot ploidy integer representing the ploidy of the data. There are
#'   ploidy * the number of samples rows in 'loci'.
#' @slot strata a factor or vector that can be coerced as to a factor as long 
#'   as the number of samples representing the current stratification scheme.
#' @slot schemes a data.frame with stratification schemes in each column.
#'   Sample names are in the rownames and must match the first part of the
#'   sample names (rownames) of the 'loci' slot. Each column is a factor.
#' @slot description a label for the object (optional).
#' @slot other a slot to carry other related information - unused in package
#'   analyses (optional).
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \code{\link{df2gtypes}}, \code{\link{sequence2gtypes}},
#'   \code{\link{genind2gtypes}}, \code{\link{gtypes.accessors}},
#'   \code{\link{initialize.gtypes}}
#' @examples
#'
#' #--- create a diploid (microsatellite) gtypes object
#' data(dolph.msats)
#' data(dolph.strata)
#' strata.schemes <- dolph.strata[, c("broad", "fine")]
#' rownames(strata.schemes) <- dolph.strata$id
#' msats.g <- new("gtypes", gen.data = dolph.msats[, -1], ploidy = 2,
#'                ind.names = dolph.msats[, 1], schemes = strata.schemes)
#' msats.g
#'
#' #--- create a haploid sequence (mtDNA) gtypes object and label haplotypes
#' data(dolph.seqs)
#' dloop.haps <- cbind(dLoop = dolph.strata$id)
#' rownames(dloop.haps) <- dolph.strata$id
#' dloop.g <- new("gtypes", gen.data = dloop.haps, ploidy = 1, 
#'                schemes = strata.schemes, sequences = dolph.seqs, 
#'                strata = "fine")
#' dloop.g
#' dloop.g <- labelHaplotypes(dloop.g, "Hap.")$gtypes
#' dloop.g
#'
#' @aliases gtypes
#' @import adegenet ape apex
#' @importFrom methods setClass
#' @export
#' 
setClass(
  Class = "gtypes",
  slots = c(loci = "data.frameOrNULL", sequences = "dnaSequences",
            ploidy = "integer", strata = "factorOrNULL",
            schemes = "data.frameOrNULL", description = "charOrNULL", 
            other = "ANY"
  ),
  prototype = prototype(loci = NULL, sequences = NULL, ploidy = 0L, 
                        strata = NULL, schemes = NULL, 
                        description = NULL, other = NULL
  ),
  validity = function(object) {
    # check that all columns in loci are factors
    loci.are.factors <- sapply(object@loci, is, class2 = "factor")
    if(!all(loci.are.factors)) {
      cat("all columns in the 'loci' slot are not factors\n")
      return(FALSE)
    }
    
    # check sequences
    if(!is.null(object@sequences)) {  
      # check that length of sequences equals number of columns in loci
      dna <- getSequences(sequences(object), simplify = FALSE)
      num.seqs <- length(dna)
      if(num.seqs > 0 & num.seqs != ncol(object@loci)) {
        cat("the number of sets of sequences is not equal to the number of loci\n")
        return(FALSE)
      }
      
      # check that locus names are the same in the @loci colnames and names of 
      #  @sequences
      if(!identical(colnames(object@loci), getLocusNames(object@sequences))) {
        cat("the names of the sets of sequences is not the same as the loci\n")
        return(FALSE)
      }
      
      # check that sequence haplotype labels can be found
      locus.good <- sapply(colnames(object@loci), function(x) {
        haps <- unique(as.character(object@loci[, x]))
        seqs <- rownames(as.matrix(dna[[x]]))
        all(na.omit(haps) %in% seqs)
      })
      if(!all(locus.good)) {
        bad.loci <- paste(colnames(object@loci)[!locus.good], collapse = ", ")
        cat("haplotypes are missing in", bad.loci, "\n")
        return(FALSE)
      }
    }
    
    # check that ploidy is compatible with loci
    if(nrow(object@loci) %% object@ploidy != 0) {
      cat("number of alleles is not an even multiple of 'ploidy'\n")
      return(FALSE)
    }
    if(object@ploidy != 1 & !is.null(object@sequences)) {
      cat("sequences can't be present unless object is haploid (ploidy = 1)\n")
      return(FALSE)
    }
    
    # check that strata is same length as number of individuals
    if(!is.null(object@strata)) {
      if((nrow(object@loci) / object@ploidy) != length(object@strata)) {
        cat("length of 'strata' slot is not equal to number of individuals\n")
        return(FALSE)
      }
    }
    
    # check that at least some individuals are in strata schemes
    if(!is.null(object@schemes)) {
      if(length(intersect(indNames(object), rownames(object@schemes))) == 0) {
        cat("no sample ids from the 'loci' slot are in stratification schemes\n")
        return(FALSE)
      }
    }
    
    TRUE
  }
)