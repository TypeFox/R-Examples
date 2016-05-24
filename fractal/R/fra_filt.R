######################################################
## FRACTAL nonlinear denoising functions
##
##   medianFilter
##   localProjection
##
######################################################

###
# medianFilter
###

"medianFilter" <- function(x, order=2)
{
  if (is(x,"signalSeries"))
    x <- x@data
  checkVectorType(x,"numeric")
  checkScalarType(order,"integer")

  as.vector(itCall( "RS_fractal_filter_median",
    x, order))
    #COPY=rep(FALSE,2),
    #CLASSES=c("matrix", "integer"),
    #PACKAGE="ifultools"))
}

###
# localProjection
###

"localProjection" <- function(x, dimension=3, tlag=timeLag(x), n.neighbor=dimension + 1,
  max.distance=2*stdev(x), metric=Inf, noise.dimension=1, corr.curve=TRUE)
{
  # check input argument types and lengths
  checkVectorType(x,"numeric")
  checkScalarType(dimension,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(n.neighbor,"integer")
  checkScalarType(metric,"numeric")
  checkScalarType(noise.dimension,"integer")
  checkScalarType(corr.curve,"logical")
  checkScalarType(max.distance,"numeric")

  # check input argument values
  if (dimension < 3)
    stop("Input dimension must be at least three")
  if (tlag < 1)
    stop("tlag must be positive")
  nembed <- length(x) - (dimension - 1) * tlag
  if (nembed < 2)
    stop("Not enough embedding points in the phase space")
  if (n.neighbor <= dimension)
    stop("Minimum number of neighbors must be greater than the embedding dimension")
  if (n.neighbor > nembed)
    stop("Input n.neighbor too large, not enough points in the tlag embedding")
  if (max.distance <= 0.0)
    stop("Input max.distance must be positive")
  if (noise.dimension < 1)
    stop("Input noise.dimension must be positive")
  if (noise.dimension >= dimension)
    stop("Input noise.dimension must be less than the embedding dimension")

  as.vector(itCall("RS_fractal_filter_nonlinear_local_projection",
    as.matrix(x),
    as.integer(dimension),
    as.integer(tlag),
    as.integer(n.neighbor),
    max.distance,
    mutilsDistanceMetric(metric),
    as.integer(noise.dimension),
    corr.curve))
    #COPY=rep(FALSE,8),
    #CLASSES=c("matrix",rep("integer",3),"numeric",rep("integer",2),"logical"),
    #PACKAGE="ifultools"))
 }
