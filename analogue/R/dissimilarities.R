###########################################################################
##                                                                       ##
## dissimilarities - Extracts dissimilarities from fitted models         ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object on which method dispatch applied (Only 'mat')  ##
## which         - For class 'analog', which dissimilarities to return   ##
##                                                                       ##
###########################################################################
dissim <- function(object, ...)
  UseMethod("dissimilarities")

dissimilarities <- function(object, ...)
  UseMethod("dissimilarities")

dissimilarities.analog <- function(object, which = c("train", "analogs"), ...)
  {
    which <- match.arg(which)
    if(which == "train") ##as.vector(as.dist(object$train))
      retval <- object$train[lower.tri(object$train)]
    else
      retval <- as.vector(object$analogs)
    class(retval) <- "dissimilarities"
    retval
  }

`dissimilarities.mat` <- function(object, ...) {
    retval <- as.vector(object$analogs)
    class(retval) <- "dissimilarities"
    retval
}
