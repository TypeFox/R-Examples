#' @title Convert a data.frame to gtypes
#' @description Load allelic data from a data.frame or matrix into a 
#'   \linkS4class{gtypes} object. 
#' 
#' @param x a matrix or data.frame of genetic data.
#' @param ploidy number of number of columns in \code{x} storing alleles
#'   at each locus.
#' @param id.col column name or number where individual sample ids are stored.
#'   If \code{NULL} then rownames are used. If there are no rownames, then 
#'   samples are labelled with consecutive numbers.
#' @param strata.col column name or number where stratification is stored. If 
#'   \code{NULL} then all samples are in one (default) stratum.
#' @param loc.col column number of first allele of first locus.
#' @param sequences a list, matrix, \code{\link{DNAbin}}, or 
#'   \linkS4class{multidna} object containing sequences. 
#' @param schemes an optional data.frame of stratification schemes.
#' @param description a label for the object (optional).
#' @param other a slot to carry other related information - unused in package
#'   analyses (optional).
#' 
#' @details
#' The genetic data in \code{x} starting at \code{loc.col} should be 
#' formatted such that every consecutive \code{ploidy} columns represent 
#' alleles of one locus. Locus names are taken from the column names in 
#' \code{x} and should be formatted with the same root locus name, with 
#' unique suffixes representing allels (e.g., for Locus1234: Locus1234.1 
#' and Locus1234.2, or Locus1234_A and Locus1234_B). \cr\cr
#' If sequences are provided in \code{sequences}, then they should be named 
#' and match haplotype labels in \code{loc.col} of \code{x}. If multiple 
#' genes are given as a \linkS4class{multidna}, then they should have the 
#' same names as column names in \code{x} from \code{loc.col} to the end.
#' 
#' @return a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{sequence2gtypes}, \link{gtypes2df},
#'   \link{gtypes2genind}, \link{gtypes2loci}
#' 
#' @examples
#' #--- create a diploid (microsatellite) gtypes object
#' data(dolph.msats)
#' ms.g <- df2gtypes(dolph.msats, ploidy = 2, strata.col = NULL, loc.col = 2)
#' ms.g
#' 
#' #' #--- create a haploid sequence (mtDNA) gtypes object
#' data(dolph.strata)
#' data(dolph.haps)
#' 
#' seq.df <- dolph.strata[ c("id", "broad", "dLoop")]
#' dl.g <- df2gtypes(seq.df, ploidy = 1, sequences = dolph.haps)
#' dl.g
#' 
#' @importFrom methods new
#' @export
#' 
df2gtypes <- function(x, ploidy, id.col = 1, strata.col = 2, loc.col = 3, 
                      sequences = NULL, schemes = NULL, description = NULL, 
                      other = NULL) {
  # check x
  if(!(is.matrix(x) | is.data.frame(x))) {
    stop("'x' must be a matrix or data.frame")
  }
  
  # check that id.col, strata.col, and loc.col are numbers and loc.col is max
  if(!is.null(id.col)) id.col <- as.numeric(id.col)
  if(!is.null(strata.col)) strata.col <- as.numeric(strata.col)
  loc.col <- as.numeric(loc.col)
  if(loc.col < max(id.col, strata.col, loc.col)) {
    stop("'loc.col' must be greater than 'id.col' and 'strata.col'")
  }
  
  # extract id, strata, and locus information
  rownames(x) 
  ind.names <- if(is.null(id.col)) {
    if(is.null(rownames(x))) {
      1:nrow(x) 
    } else {
      rownames(x)
    }
  } else x[, id.col]
  ind.names <- as.character(ind.names)
  strata <- if(is.null(strata.col)) NULL else x[, strata.col]
  gen.data <- x[, loc.col:ncol(x), drop = FALSE]
  
  # return new gtypes object
  new("gtypes", gen.data = gen.data, ploidy = ploidy, ind.names = ind.names,
      strata = strata, schemes = schemes, sequences = sequences, 
      description = description, other = other
  )
}