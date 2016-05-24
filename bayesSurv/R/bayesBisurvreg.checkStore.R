###################################################
#### AUTHOR:     Arnost Komarek                ####
####             (2005)                        ####
####                                           ####
#### FILE:       bayesBisurvreg.checkStore.R   ####
####                                           ####
#### FUNCTIONS:  bayesBisurvreg.checkStore     ####
###################################################

### ======================================
### bayesBisurvreg.checkStore
### ======================================
## 21/12/2004: start working on it
##
bayesBisurvreg.checkStore <- function(store)
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

  ia2 <- match("a2", instore, nomatch=NA)
  if(is.na(ia2)) store$a2 <- FALSE
  if (!is.logical(store$a2)) store$a2 <- FALSE

  iy2 <- match("y2", instore, nomatch=NA)
  if(is.na(iy2)) store$y2 <- FALSE
  if (!is.logical(store$y2)) store$y2 <- FALSE

  ir2 <- match("r2", instore, nomatch=NA)
  if(is.na(ir2)) store$r2 <- FALSE
  if (!is.logical(store$r2)) store$r2 <- FALSE

  return(store)  
}  
