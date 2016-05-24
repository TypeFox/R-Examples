#' @title Trim N's From Sequences
#' @description Removes N's from beginning and end of sequences.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' 
#' @return sequences with beginning and trailing N's removed.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
trimNs <- function(x) {
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  
  dna <- as.character(as.list(x))
  result <- lapply(1:length(dna), function(i) {
    seq.vec <- paste(dna[[i]], collapse = "")
    start <- gregexpr("^[n]+", seq.vec)[[1]]
    end <- gregexpr("[n]+$", seq.vec)[[1]]
    start <- ifelse(start == -1, 1, attr(start, "match.length") + 1)
    end <- ifelse(end == -1, nchar(seq.vec), end - 1)
    dna[[i]][start:end]
  })
  
  result <- as.DNAbin(result)
  names(result) <- names(x)
  result
}