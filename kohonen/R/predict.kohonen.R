"predict.kohonen" <- function(object, newdata, trainX, trainY,
                              unit.predictions = NULL, threshold = 0,
                              whatmap = NULL, weights = 1, ...)
{
  mapping <- NULL
  
  if (missing(newdata)) {
    if (!is.null(object$data)) {
      newdata <- object$data
      mapping <- object$unit.classif # perhaps NULL
    } else {
      stop("No data given with which to predict")
    }
  }

  ## Find the mapping for the new data
  if (is.null(mapping))
    mapping <- map.kohonen(object, newdata, whatmap, weights)$unit.classif
  
  ## Find the value for each unit. For unsupervised maps, we should
  ## have trainX and trainY values. Argument unit.predictions will
  ## override the predictions per unit.
  if (is.null(unit.predictions)) {
    if (object$method %in% c("xyf", "bdk")) {
      unit.predictions <- list(object$codes$Y)
      contin <- object$contin
      factorY <- !contin
      if (factorY)
          factorY.levels <- list(colnames(object$codes$Y))
    } else {
      ## If whatmap is given, and not all layers used in training the map are
      ## included, the excluded layers are interpreted as being the ones
      ## for which a prediction is wanted. Unit.predictions are already
      ## available.
      if (object$method == "supersom" && !is.null(whatmap)) {
        whatmap <- check.whatmap(object, whatmap)
        
        trained.layers <- object$whatmap[!(object$whatmap %in% whatmap)]
        if (length(trained.layers) > 0) {
          unit.predictions <- object$codes[trained.layers]
          contin <- object$contin[trained.layers]
          factorY <- !contin
          if (any(factorY))
              factorY.levels <-
                  lapply(object$codes[trained.layers[factorY]], colnames)
        }
      }

      ## If unit.predictions are not available, we have to find
      ## them... For SOMs and supersoms.
      if (is.null(unit.predictions)) {
        if (missing(trainY))
          stop("For unsupervised forms of mapping, trainY is required")
        if (!is.list(trainY))
          trainY <- list(trainY)
        factorY <- sapply(trainY, is.factor)
        factorY.levels <- lapply(trainY[factorY], levels)
        trainY[factorY] <- lapply(trainY[factorY], classvec2classmat)

        trainY[sapply(trainY, is.vector)] <-
          lapply(trainY[sapply(trainY, is.vector)], matrix, ncol = 1)

        contin <- sapply(trainY,
                         function(x) any(abs(rowSums(x) - 1) > 1e-8))
        nY <- sapply(trainY, ncol)
        
        trainingMapping <- NULL
        if (missing(trainX) & !is.null(object$data)) {
          trainX <- object$data
          trainingMapping <- object$unit.classif
        }

        nX <- ifelse(is.list(trainX),
                     nrow(trainX[[1]]),
                     nrow(trainX))
        if (nX != nrow(trainY[[1]]))
          stop("Unequal number of rows in trainX and trainY")
        
        ## Find mapping for training data
        if (is.null(trainingMapping))
          trainingMapping <- map.kohonen(object, trainX)$unit.classif
        
        nhbrdist <- unit.distances(object$grid, object$toroidal)

        ## Find unit.predictions for training data; loop over list elements
        unit.predictions <- vector(length(nY), mode = "list")
        names(unit.predictions) <- names(trainY)
        
        for (ii in seq(along = trainY)) {
          unit.predictions[[ii]] <- matrix(NA, nrow(object$grid$pts), nY[ii])
          huhn <- aggregate(trainY[[ii]], by = list(cl = trainingMapping),
                            mean)
          unit.predictions[[ii]][huhn[,1],] <- as.matrix(huhn[,-1])
          
          ## Prediction for empty units
          nas <- which(apply(unit.predictions[[ii]],
                             1,
                             function(x) all(is.na(x))))
          for (i in seq(along = nas)) {
            unit.predictions[[ii]][nas[i], ] <-
              colMeans(unit.predictions[[ii]][nhbrdist[nas[i],] == 1, ,
                                              drop=FALSE], na.rm = TRUE)
          }
          
          colnames(unit.predictions[[ii]]) <- colnames(trainY[[ii]])
        }
      }
    }
  }

  prediction <- lapply(unit.predictions, function(x) x[mapping,,drop = FALSE])
  if (any(!contin))
      prediction[!contin] <- lapply(prediction[!contin],
                                    classmat2classvec, threshold = threshold)

  for (fY in which(factorY))
      prediction[[fY]] <- factor(prediction[[fY]],
                                 levels = factorY.levels[[fY]])

  if (length(prediction) == 1) {
    unit.predictions <- unit.predictions[[1]]
    prediction <- prediction[[1]]
  }
  
  list(prediction = prediction, unit.classif = mapping,
       unit.predictions = unit.predictions)
}
