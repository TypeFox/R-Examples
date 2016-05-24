# file:    sortFrame.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# sortFrame() sorts the cases of a data frame by one or more variables. There's a bit of 
# work to do to make the sorting more intuitive, but it's fairly functional now.
sortFrame <- function(x,..., alphabetical = TRUE){
  
  if( !is(x,"data.frame") ) {
    stop( '"x" must be a data frame')
  }
  if( !is(alphabetical,"logical") | length(alphabetical) !=1 ) {
    stop( '"alphabetical" must be a single logical value')
  }
  
  dots <- as.list(substitute(list(...)))[-1L] # list of quoted sort terms
  if( length(dots) == 0 ){ return(x) } # do nothing if null arguments
  rel.vars <- unlist(lapply(dots,all.vars)) # which variables are referred to
  y <- lapply(x[rel.vars], xtfrm) # numeric frame that sorts equivalently
  if( alphabetical == TRUE ) { # case conversion if necessary...   
    char.vars <- unlist(lapply(x,is,"character")) # find character vars
    char.vars <- names(which(char.vars)) # relevant variable names
    n <- length(y[[1]]) # number of cases
    for( v in char.vars ) { 
      z <- xtfrm(tolower(x[[v]])) # sorts equivalently to lower case version 
      y[[v]] <- z + y[[v]]/(n+1) # original only breaks ties
    }
  }
  
  ord <- with(y, do.call(order, dots)) # the sort order
  return( x[ord,] ) # sort and return
  
}