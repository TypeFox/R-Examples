
##  Function stripNAs
##
## Remove NA (missing data) entries from a vector.
##
## Removes NA entries from a vector.  R already contains a function,
## \code{\link{na.omit}} for this purpose, but for some reason the
## result from this function breaks some of the PPE code.
##
## This function is a utility function provided primarily for the
## benefit of other functions within the PPE package.
##
## @title Remove NA's from a vector
## @param vv a vector.
## @return The input vector, excluding any NA entries.
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## vv <- c(1:3,NA,5:7,NA,NA)
## print(vv)
## stripNAs(vv)
##
stripNAs <- function(vv) { return(vv[!is.na(vv)]) }
