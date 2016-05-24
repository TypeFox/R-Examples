################################################################################
##
## $Id: nearest.multiple.R 346 2006-10-01 05:08:55Z enos $
##
## Find the multiple of y nearest to x.
##
################################################################################

.nearest.multiple <- function(x, y){

  ## Currently only defined for y != 0.
  
  stopifnot(all(y != 0, na.rm = TRUE))

  y <- abs(y)
  y <- ifelse(x < 0, y * -1, y)

  ## Take modulus wrt rount.lot, then round up or down to the nearest
  ## lot.
  
  x.mod.mult <- x %% y

  x <- ifelse(abs(x.mod.mult) >= abs(y/2),
              x + y - x.mod.mult,
              x - x.mod.mult)
  x
}
