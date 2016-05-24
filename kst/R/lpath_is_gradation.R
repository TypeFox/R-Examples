### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### lpath_is_gradation.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2009-04-20: created
###

lpath_is_gradation <- function(x) {

   ### check for gradation
   FUN <- function(y) {
      sd <- diff(unlist(lapply(y, function(z) length(z))), lag=1)
      all(sd <= 1)
   }
   lapply(x, FUN)

}
