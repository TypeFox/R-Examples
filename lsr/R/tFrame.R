# file:    tFrame.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 29 June 2013

# tFrame() transposes a data frame. This isn't usually a very sensible thing to do, 
# except in those cases where the data frame would actually make sense as a matrix,
# in which case you could coerce to a matrx and then transpose using t(). However,
# there are some cases where novices (who haven't yet grasped coercion) might want 
# to transpose a data frame, and it can be handy to use the tFrame() function for 
# teaching purposes rather than stop to teach coercion & data types on the spot.
tFrame <- function(x) {
  if (!is(x,"data.frame")) {
    stop("'tFrame' is intended to apply to data frames only")
  }
  x <- t(x) # coerce to matrix and transpose
  x <- as.data.frame(x)
  return( x )
}