# file:    permuteLevels.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# permuteLevels() allows the user to reorder the levels of a factor according 
# to some arbitrary permutation. It's much like relevel() but more general. 
permuteLevels <- function(x,perm,ordered = is.ordered(x),invert=FALSE) {
  
  if( !is(x,"factor") ) stop( '"x" must be a factor')
  if( !is(perm,"numeric") ) stop('"perm" must be numeric') 
  if( length(perm) != nlevels(x)) stop('length of "perm" must equal the number of levels of "x"')
  if( any( (sort(perm) - 1:nlevels(x)) != 0)) stop( '"perm" is not a valid permutation')  
  if( !is(ordered,"logical") | length(ordered) !=1 ) {
    stop( '"ordered" must be a single logical value')
  }
  if( !is(invert,"logical") | length(invert) !=1 ) {
    stop( '"invert" must be a single logical value')
  }
    
  if(invert){ perm <- order(perm) }
  y <- factor( x = order(perm)[as.numeric(x)],
               levels = seq_along(levels(x)),
               labels = levels(x)[perm],
               ordered = ordered ) 
  return(y)
}