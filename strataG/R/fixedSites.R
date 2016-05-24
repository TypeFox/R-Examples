#' @title Fixed Sites
#' @description Identify fixed sites among sequences.
#' 
#' @param x a \code{\link{gtypes}} object with sequences, a list of sequences, or a consensus sequence. 
#'   Sequences must be aligned.
#' @param bases a character vector of valid bases to consider.
#' 
#' @return a vector of fixed sites. Element names are site positions in the
#'   original sequence.
#' 
#' @author Eric Archer <eric.archer@@noaa.gov>
#' 
#' @examples
#' data(dolph.haps)
#' 
#' fixedSites(dolph.haps)
#' 
#' @seealso \code{\link{variableSites}}
#' 
#' @export
#' 
fixedSites <- function(x, bases = c("a", "c", "g", "t", "-")) {
  x <- as.multidna(x)
  bases <- tolower(bases)
  
  result <- lapply(getSequences(x, simplify = FALSE), function(dna) {
    dna <- as.character(as.matrix(dna))
    is.fixed <- apply(dna, 2, function(site) {
      site <- site[site %in% bases]
      length(unique(site)) == 1
    })
    sites <- which(is.fixed)
    fixed.bp <- sapply(sites, function(i) {
      bps <- dna[, i]
      unique(bps[bps %in% bases])[1]
    })
    names(fixed.bp) <- sites
    fixed.bp
  })
  
  if(length(result) == 1) {
    result[[1]]
  } else {
    names(result) <- getSequenceNames(x)
    result
  }
}