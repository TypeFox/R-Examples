# file:    maxFreq.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 18 February 2015

# maxFreq() returns the frequency of the sample mode.
maxFreq <- function(x, na.rm = TRUE) {
  
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
  max.freq <- max(freq, na.freq)                      # modal frequency    
  return( max.freq )                                 
  
}
