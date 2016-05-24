### predict.wccsom is a function to predict properties for SOMs and
### XYFs; i.e. the average of the property under investigation for
### all objects mapped on that particular unit.
### Arguments trainX, trainY and unit.predictions are optional: when
### no training data have 
### been stored in the network, or when (in the supervised case) a
### complicated Y-structure is present.

predict.wccsom <- function(object, newdata, 
                           trainX, trainY, unit.predictions, ...)
{
  mapping <- NULL
  
  if (missing(newdata)) {
    if (!is.null(object$data)) {
      newdata <- object$data
      mapping <- object$unit.classif
    } else {
      stop("No data given with which to predict")
    }
  }

  ## unit.predictions are the properties of interest for the units in
  ## the map
  if (missing(unit.predictions)) {
    if (is.null(object$predict.type)) {
      ## Prediction for SOMs
      ## First calculate averages of properties for every unit. Empty
      ## units are estimated from neighbouring values by interpolation.
      if (missing(trainX) & !is.null(object$data))
        trainX <- object$data
      
      if (!missing(trainY)) {
        if (is.vector(trainY)) trainY <- matrix(trainY, ncol=1)
        nY <- ncol(trainY)
      } else {
        stop("For unsupervised maps, argument trainY is required")
      }
      
      if (is.null(object$unit.classif))
        object$unit.classif <- wccassign(object, trainX)

      unit.predictions <- matrix(NA, nrow(object$codes), nY)
      huhn <- aggregate(trainY,
                        by=list(cl = object$unit.classif),
                        mean)
      
      ## the next line should be simpler in a way...
      unit.predictions[sort(as.numeric(levels(huhn[,1]))),] <-
        as.matrix(huhn[,-1])
      
      nas <- which(apply(unit.predictions, 1, function(x) all(is.na(x))))
      nhbrdist <- unit.distances(object$grid, object$toroidal)
      for (i in seq(along=nas)) {
        unit.predictions[nas[i],] <-
          mean(unit.predictions[nhbrdist[i,] == 1,],
               na.rm=TRUE)
      }

      colnames(unit.predictions) <- colnames(trainY)
    } else {
      unit.predictions <- object$codeYs
    }
  }

  if (is.vector(unit.predictions))
    unit.predictions <- matrix(unit.predictions, ncol=1)

  ## After having determined the unit.predictions, we only have to map
  ## the new data and then we can return the unit.predictions
  if (is.null(mapping))
    mapping <- wccassign(object, newdata)$classif

  list(unit.predictions = unit.predictions,
       predictions = unit.predictions[mapping,])
}
