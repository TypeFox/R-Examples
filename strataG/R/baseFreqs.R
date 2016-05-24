#' @title Base Frequencies
#' @description Calculate nucleotide base frequencies along a sequence.
#' 
#' @param x a \linkS4class{gtypes} object with aligned sequences or a list of 
#'   aligned DNA sequences.
#' @param bases character vector of bases. Must contain valid IUPAC codes. 
#'   If \code{NULL}, will return summary of frequencies of observed bases.
#' @param ignore a character vector of bases to ignore when 
#'   calculating frequencies.
#' 
#' @return For each gene, a list containing:
#' \tabular{ll}{
#'   \code{site.freqs} \tab a matrix of base frequencies at each site.\cr
#'   \code{base.freqs} \tab a vector of overall base frequencies.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
baseFreqs <- function(x, bases = NULL, ignore = c("n", "x", "-", ".")) {
  
  bases <- if(is.null(bases)) {
    rownames(iupac.mat)
  } else {
    tolower(as.character(bases))
  }
  ignore <- tolower(ignore)
  
  x <- getSequences(as.multidna(x), simplify = FALSE)
  result <- lapply(x, function(dna) {
    seq.mat <- as.character(as.matrix(dna))
    site.freqs <- apply(seq.mat, 2, function(site) {
      site <- site[!site %in% ignore]
      table(factor(site, levels = bases))
    })
    colnames(site.freqs) <- 1:ncol(site.freqs)
    base.freqs <- table(factor(as.vector(seq.mat), levels = bases))
    base.freqs <- base.freqs / sum(base.freqs)
    list(site.freqs = site.freqs, base.freqs = base.freqs)
  })
  
  if(length(result) == 1) {
    result[[1]]
  } else {
    names(result) <- names(x)
    result
  }
}