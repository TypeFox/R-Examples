#' @title Compiles and Formats Allele Frequencies
#' @description Format allele frequencies to for a set of ids and loci 
#' 
#' @param x a matrix or data.frame where first column is sample id and 
#'   second colum is locus name.
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return data.frame of original samples, loci, and formatted alleles and 
#'   frequencies.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' x <- cbind(
#'  id = sample(indNames(msats.g), 10, rep = TRUE),
#'  locus = sample(locNames(msats.g), 10, rep = TRUE)
#' )
#' alleleFreqFormat(x, msats.g)
#' 
#' @export
#' 
alleleFreqFormat <- function(x, g) {
  if(!(is.data.frame(x) | is.matrix(x))) {
    stop("'x' must be a data.frame or matrix")
  }
  if(ncol(x) != 2) stop("'x' must have two columns")
  x <- as.matrix(x)
  
  freqs <- alleleFreqs(g)
  fmtd <- rep(as.character(NA), nrow(x))
  for(i in 1:length(fmtd)) {
    id <- x[i, 1]
    locus <- x[i, 2]
    # skip (leave as NA) if either id or locus can't be found
    if(!(id %in% indNames(g) | locus %in% locNames(g))) next
    # get genotype of this id at this locus
    gt <- unlist(loci(g, id, locus))
    # if the genotype is NA skip and leave format as NA
    if(any(is.na(gt))) next
    # get frequency and round
    f <- round(freqs[[locus]][gt, "prop"], 3)
    f <- paste(gt, " (", format(f, nsmall = 3), ")", sep = "")
    # return a single frequency if homozygote
    fmtd[i] <- if(length(unique(gt)) == 1) f[1] else paste(f, collapse = " / ")
  }
  cbind(x, allele.freqs = fmtd)
}