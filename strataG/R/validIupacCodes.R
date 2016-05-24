#' @title Valid IUPAC Codes
#' @description Get all possible valid IUPAC codes.
#' 
#' @param bases character vector of nucleotides or IUPAC codes to be checked.
#' 
#' @return character vector of valid IUPAC codes for \code{bases}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' validIupacCodes(c("c", "t", "c", "c"))
#' 
#' validIupacCodes(c("c", "y", "c", "c"))
#' 
#' validIupacCodes(c("a", "g", "t", "a"))
#' 
#' @export
#' 
validIupacCodes <- function(bases) {
  bases <- tolower(bases)
  base.rows <- which(rownames(iupac.mat) %in% bases)
  if(length(base.rows) == 0) stop("No valid IUPAC codes in 'bases'")
  valid.codes <- sapply(colnames(iupac.mat), 
                        function(code) all(iupac.mat[base.rows, code]))
  colnames(iupac.mat)[valid.codes]
}