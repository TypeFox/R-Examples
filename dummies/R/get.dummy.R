# -----------------------------------------------------------------------------
# FUNCTIONS: get.dummy get.dummies
#  Retrieve the dummy variables for a data.frame
# -----------------------------------------------------------------------------

get.dummy <- function(data,name=NULL) {
    
  if( ! is.null(name) ) {
    dat <- data[ , which.dummy(data, name) ]  
  } else {
    dat <- data[ , which.dummy(data) ] 
  }

  return(dat)
}
