checkSetup = function(object, quiet = FALSE) {
	if (!quiet)
		cat("Checking object ... ")
# check observations
	if (is.null(object$observations))
		stop("no observations provided")
	if (!is(object$observations, "Spatial") && 
      "data" %in% names(getSlots(class(object$observations))))
		stop("observations not of class Spatial*DataFrame")
	if (dim(coordinates(object$observations))[1] < 20)
		stop("Less than 20 observations provided, not able to perform interpolation")
  if (sum(duplicated(coordinates(object$observations)))>0)
    stop("Some of the observation locations are duplicated, not able to perform interpolation")


# Check formulaString
	if (!"formulaString" %in% names(object)) 
	  stop("object does not contain a formulaString")
  depVar=as.character(object$formulaString[[2]])
	if (is.null(object$observations[[depVar]]))
		stop(paste("observations has no attribute with name",depVar))

# Check that simulations are not requested for wrong method
  if ("nsim" %in% names(object$outputWhat) && "copula" %in% class(object))
    stop("The copula method is not able to return simulations, suggesting to change method to transGaussian")

# check PredictionLocations
	if (is.null(object$predictionLocations))
		stop("no predictionsLocations provided")
	if (!is(object$predictionLocations, "Spatial"))
		stop("predictionLocations should be of class extending Spatial")

# Check if targetCRS was set
  if (!is.null(object$targetCRS)) {
# check if observations CRS was set
	  if (is.na(proj4string(object$observations)))
		  stop("can't reproject observations when its CRS is not set")
# check if PredictionLocations CRS was set:
	  if (is.na(proj4string(object$predictionLocations)))
		  stop("can't reproject predictionLocations when CRS is not set")
  }
# check that the biases to be added is among the ones removed:
  if (!is.na(object$params$addBias) & 
                sum(object$params$addBias %in% object$params$removeBias) 
                                        < length(object$params$addBias))
    stop("Cannot add biases that have not been removed")

	if (dim(coordinates(object$predictionLocations))[1] > 2  &&
     commonArea(object$observations,object$predictionLocations)[2] < 0.01)
	  warning("Boundary boxes of predictionLocations and observations show small overlap, check projections if predictionLocations are not supposed to be concentrated")
	if (!quiet)
		cat("OK\n")
	invisible(TRUE);
}



