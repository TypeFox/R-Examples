#' @title Sequence Summaries
#' @description Summaries for each sequence.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' 
#' @return a matrix listing the start and end positions of each sequence 
#'   (excluding beginning and trailing N's), the length, the number of N's, 
#'   and the number of indels.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' library(apex)
#' data(woodmouse)
#' 
#' summarizeSeqs(woodmouse)
#' 
#' @export
#' 
summarizeSeqs <- function(x) {
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  
  x <- sapply(
    x, function(dna) unlist(as.character(as.list(dna))), simplify = FALSE
  )

  t(sapply(x, function(this.seq) {
    seq.rle <- rle(this.seq)
    start <- ifelse(seq.rle$values[1] == "-", seq.rle$lengths[1] + 1, 1)
    start <- unname(start)
    n <- length(seq.rle$values)
    sum.n <- sum(seq.rle$lengths)
    end <- ifelse(seq.rle$values[n] == "-", sum.n - seq.rle$lengths[n], sum.n)
    end <- unname(end)
    num.ns <- sum(this.seq == "n")
    num.indels <- sum(this.seq[start:end] == "-")
    c(start = start, end = end, length = end - start + 1, 
      num.ns = num.ns, num.indels = num.indels
    )
  }))
}
