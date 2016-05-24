#' @name parseLocusNames
#' @title Parse Locus Names
#' @description Parse locus names from a matrix where alleles are in
#'   consecutive columns.
#'
#' @param locus.names a vector of column names, where each locus must be
#'   named with the same roots. For example, diploid locus 'ABCD' would have
#'   two columns named something like 'ABCD.1' and 'ABCD.2', or
#'   'ABCD_A' and 'ABCD_B'.\cr
#' @param ploidy integer representing the ploidy of the data.
#'
#' @return A vector of the root locus names.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
.parseLocusNames <- function(locus.names, ploidy) {
  if(ploidy == 1) return(locus.names)
  loc.i <- matrix(1:length(locus.names), nrow = ploidy)
  apply(loc.i, 2, function(i) {
    this.loc <- locus.names[i]
    max.length <- max(sapply(this.loc, nchar))
    # find location of first difference
    ptr <- 1
    all.same <- length(unique(substr(this.loc, 1, ptr))) == 1
    while(all.same & ptr < max.length) {
      ptr <- ptr + 1
      all.same <- length(unique(substr(this.loc, 1, ptr))) == 1
    }
    ptr <- ptr - 1
    if(ptr == 0) {
      this.loc[1] # return first name if all names are different
    } else {
      # remove last same character if it is not alphanumeric
      if(!substr(this.loc[1], ptr, ptr) %in% c(LETTERS, letters, 0:9)) {
        ptr <- ptr - 1
      }
      substr(this.loc[1], 1, ptr)
    }
  })
}