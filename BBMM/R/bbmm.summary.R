bbmm.summary <-
function(x){
#
#   print method for BBMM objects
#
    cat("\nBrownian motion variance : ", x[[1]], fill=TRUE) 
    cat("Size of grid : ", length(x$x), "cells", fill=TRUE)
    cat( "Grid cell size : ", abs(x$x[1]-x$x[2]), fill=TRUE)

}
