####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       bayessurvreg1.checkStore.R     ####
####                                            ####
#### FUNCTIONS:  bayessurvreg1.checkStore       ####
####################################################

### ======================================
### bayessurvreg1.checkStore
### ======================================
bayessurvreg1.checkStore <- function(store)
{
  if(length(store) == 0) instore <- "arnost"
  else                   instore <- names(store)

  iy <- match("y", instore, nomatch=NA)
  if(is.na(iy)) store$y <- FALSE
  if (!is.logical(store$y)) store$y <- FALSE

  ir <- match("r", instore, nomatch=NA)
  if(is.na(ir)) store$r <- FALSE
  if (!is.logical(store$r)) store$r <- FALSE

  ib <- match("b", instore, nomatch=NA)
  if(is.na(ib)) store$b <- FALSE
  if (!is.logical(store$b)) store$b <- FALSE

  iu <- match("u", instore, nomatch=NA)
  if(is.na(iu)) store$u <- FALSE
  if (!is.logical(store$u)) store$u <- FALSE

  iMHb <- match("MHb", instore, nomatch=NA)
  if(is.na(iMHb)) store$MHb <- FALSE
  if (!is.logical(store$MHb)) store$MHb <- FALSE

  irr <- match("regresres", instore, nomatch=NA)
  if(is.na(irr)) store$regresres <- FALSE
  if (!is.logical(store$regresres)) store$regresres <- FALSE
  
  return(store)  
}  
