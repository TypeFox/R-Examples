#
# replace NA values with "BLANK" for factor variables
#
# Author: Justin Hemann

#
# define generic function
#
CleanNaFromFactor <- function(x, ...) {
  UseMethod("CleanNaFromFactor", x)
}


CleanNaFromFactor.factor = function(x, replacement="BLANK", ...) {
  naCount <- sum( is.na(x) ) # count of occurences of NA
  if (naCount > 0) {
    # add the "BLANK" factor
    x <- factor(x, levels = c(levels(x),replacement))
    # replace NAs with BLANK in this variable
    x[is.na(x)] <- replacement
  }
  # return the factor
  return(x)
}
