#' @title IUPAC Code
#' @description Calculate the correct IUPAC code for a vector of nucleotides.
#'
#' @param bases character vector containing valid nucleotides or IUPAC codes.
#' @param ignore.gaps logical. Ignore gaps at a site when creating consensus. 
#'   If true, then bases with a gap are removed before consensus is calculated. 
#'   If false and a gap is present, then the result is a gap.
#'
#' @return a character representing the correct IUPAC code.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \code{\link{validIupacCodes}}
#' 
#' @examples
#' iupacCode(c("a", "a", "g"))
#' 
#' iupacCode(c("t", "c", "g"))
#'
#' @export
#' 
iupacCode <- function(bases, ignore.gaps = FALSE) {
  bases <- as.character(bases)
  if(ignore.gaps) bases <- bases[!bases %in% c("-", ".")]
  validIupacCodes(bases)[1]
}