######################################################
## FRACTAL nearest neighbor functionality
##
##   findNeighbors
##
######################################################

###
# findNeighbors
###

"findNeighbors" <- function(x, n.neighbor=NULL, metric=1,
  max.distance=0.0, olag=0, sort.distances=TRUE)
{
  # check input argument types and lengths
  checkScalarType(metric,"numeric")
  checkScalarType(max.distance,"numeric")
  checkScalarType(olag,"integer")
  checkScalarType(sort.distances,"logical")

  # perform argument checks
  if (!is.matrix(x))
    stop("Input must be a matrix")

  # use n.neighbor if specified.
  # otherwise, set default values
  if (!is.null(n.neighbor)){

    if (n.neighbor <= 0)
      stop("n.neighbor input must be positive")

    max.distance <- 0.0
  }
  else if (max.distance <= 0.0){

    n.neighbor   <- 1
    max.distance <- 0.0
  }
  else
    n.neighbor <- 0

  checkScalarType(n.neighbor,"integer")
  if (olag < 0)
    stop("olag input must be non-negative")

  if (n.neighbor > nrow(x))
    stop("The number of nerest neighbors to find, argument n.neighbor, must be ",
          "less than the number of rows in the embedding matrix, argument x")

  # call the nearest neighbor routine
  z <- itCall( "RS_fractal_neighbor_find",
    as.matrix(x),
    as.integer(n.neighbor),
    as.numeric(max.distance),
    mutilsDistanceMetric(metric),
    as.logical(sort.distances),
    as.integer(olag))
    #COPY=rep(FALSE,6),
    #CLASSES=c("matrix", "integer", "numeric",
    #  "integer", "logical", "integer"),
    #PACKAGE="ifultools")

  # map index base 0 to index base 1
  z[[1]] <- z[[1]] + 1
  z[[2]] <- z[[2]] + 1

  names(z) <- c("original", "neighbor", "distance")

  z
}
