#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       bayesHistogram.checkStore.R         ####
####                                                 ####
#### FUNCTIONS:  bayesHistogram.checkStore           ####
#########################################################

### ======================================
### bayesHistogram.checkStore
### ======================================
bayesHistogram.checkStore <- function(store)
{
  if(length(store) == 0) instore <- "arnost"
  else                   instore <- names(store)

  ia <- match("a", instore, nomatch=NA)
  if(is.na(ia)) store$a <- FALSE
  if (!is.logical(store$a)) store$a <- FALSE

  iy <- match("y", instore, nomatch=NA)
  if(is.na(iy)) store$y <- FALSE
  if (!is.logical(store$y)) store$y <- FALSE

  ir <- match("r", instore, nomatch=NA)
  if(is.na(ir)) store$r <- FALSE
  if (!is.logical(store$r)) store$r <- FALSE

  return(store)  
}  
