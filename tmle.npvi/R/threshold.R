threshold <- function(xx, min=0, max=1-min) {
  min <- Arguments$getNumeric(min);
  max <- Arguments$getNumeric(max);
  xx <- Arguments$getNumerics(xx);
  pmin(max, pmax(min, xx))
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

