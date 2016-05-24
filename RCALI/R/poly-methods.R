# +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AB:  oct 2005, sept 2009
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# METHODS 'poly'
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# -----------------------------------------------------
# crpoly:          
# FUNCTION
# Create interactivly an object 'poly' by means of the function
# getpoly from splancs. A call to a plot is required before.
# -----------------------------------------------------
crpoly <- function() {
# Load 'splancs' if not already done: (not useful, splancs is in Depends)
#   if (!any(search() == "package:splancs"))
#     library(splancs)
 a<-getpoly()
#  Click on the requested vertices with the button 1 (left)
#  of the mouse. End  with the button 2
  a=as.poly(a)
  return(a)
} # end crpoly


# -----------------------------------------------------
# as.poly: 
# FUNCTION
# Create an object 'poly' from:
# - either 2 vectors which contains respect.
# the x- and the y-coordinates of a polygon
# - or a 2-columns matrix which contains respect.
# the x- and the y-coordinates of a polygon
# ARGUMENTS:
# -  either 2 numerical vectors of same length
# -  or a matrix with at least 2 columns, the first one contains
# the x-coordinates and the second one contains
# the y-coordinates
# ----------------------------------------
# NOTE for the user:
# To create interactivly a matrix with coordinates, by clicking
# on points of a graph, use 'getpoly' from package splancs,
# then, use this function to convert it into a 'poly' object. 
# Example:
# > library(splancs)
# > plot(x=c(1,10), y=c(1,10))
# > a<-getpoly()
#  Click on the requested vertices with the button 1 (left) 
# of the mouse. End  with the button 2
# > a=as.poly(a)
# -----------------------------------------------------
as.poly <- function(x, y=NULL) {
  if (is.null(y)) {
# A matrix
    if (ncol(x) < 2)
      stop("Argument should be a matrix with at least two columns.")
    retour <- matrix(c(x[,1], x[,2]), ncol=2)
  }

  else {
# A vector with the x and a vector with the y
    if (length(x) != length(y)) 
      stop("The x-vector and the y-vector have not the same length.")
    retour <- matrix(c(x, y), ncol=2)
  }

  if (nrow(retour) <3) 
    warning("The polygon has less than three points.")

  dimnames(retour) <- list(NULL, c("xcoord", "ycoord"))
  class(retour) <- "poly"
  return(retour)
}



# -----------------------------------------------------
# plot.poly:
# FUNCTION:
#  Plot a 'poly' object by means of the function 'polymap'
#  from splancs
# If the 'poly' has an attribute named 'couleur', it is
# plotted in this color
# ARGUMENTS:
# - x: an object of classe 'poly'
# - ...: a variable list of arguments which will be passed as it to
#  'polymap'
# -----------------------------------------------------
plot.poly <- function(x, ...) {
# We unclass the object because polymap call 'plot'
# on the object and this would call automatically 'plot.poly' 
# instead of the standard plot.
  poly <- unclass(x)

# We load 'splancs' if not already done:  (not useful, splancs is in Depends)
#   if (!any(search() == "package:splancs"))
#     library(splancs)

  if (!is.null(attr(poly, "couleur")))
    polymap(poly, col=attr(poly, "couleur"), ...)
  else
    polymap(poly, ...)
}
# -----------------------------------------------------
