# This file is intented to be called by calex_1d.R.  It defines an
# extractor function for D1.

extractor.1d <- function(D1){
return(list(x.star = D1[, 1, drop = FALSE], t.vec = D1[, 2, drop = FALSE]))
}
