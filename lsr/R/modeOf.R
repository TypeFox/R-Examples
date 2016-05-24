# file:    modeOf.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 18 February 2015

# modeOf() returns the sample mode: the value that has highest observed frequency
modeOf <- function(x, na.rm = TRUE) {
  
  if( !is.vector(x) & !is.factor(x) ) {
    stop( '"x" must be a vector or a factor')
  }
  if( !is(na.rm,"logical") | length(na.rm) !=1 ) {
    stop( '"na.rm" must be a single logical value')
  }

  na.freq <- 0                                        
  if (na.rm == FALSE) { na.freq <- sum( is.na(x) ) }  # count the NAs if needed
  x <- x[!is.na(x)]                                   # delete NAs  
  obs.val <- unique(x)                                # find unique values
  valFreq <- function(x, y){ sum(y == x) }
  freq <- unlist((lapply( obs.val, valFreq, x )))     # apply for all unique values
  max.freq <- max(freq)                               # modal frequency
  if (na.rm == FALSE & na.freq > max.freq) {
    modal.values <- NA                                # mode is NA if appropriate...
  } else {
    modal.cases <- freq == max.freq                   # otherwise find modal cases
    modal.values <- obs.val[modal.cases]              # and corresponding values
  }
  if(class(x)=="factor") modal.values <- as.character(modal.values)
  return( modal.values )
  
}