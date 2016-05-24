# file:    aad.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# aad() returns the average (mean) absolute deviation from the sample mean.
# I'd have called this function mad() but that's taken by a slightly different
# function
aad <- function(x, na.rm = FALSE) { 
  if ( !is(x,"numeric") & !is(x,"integer") ) {
    stop( '"x" must be numeric')
  }
  if( !is(na.rm,"logical") | length(na.rm) !=1 ) {
    stop( '"na.rm" must be a single logical value')
  }
  if (na.rm) { x <- x[!is.na(x)] }
  y <- mean( abs(x - mean(x)) )
  return(y)
}

