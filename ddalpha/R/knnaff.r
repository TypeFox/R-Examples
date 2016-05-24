knnaff.train <- function(data, aggregation.method = "majority", range = -1, k = -1, i = 0){

  knnaff <- knaff.create.structure(data)
  # Checks
  if (!is.character(aggregation.method) 
      || length(aggregation.method) != 1 
      || !(aggregation.method %in% c("majority", "sequent"))){
    warning("In treatment number ", i, ": Argument \"aggregation.method\" not 
            specified correctly. \"majority\" is used as a default value")
    knnaff$methodAggregation <- "majority"
  }else{
    knnaff$methodAggregation <- aggregation.method
  }
  if (!is.numeric(range) 
      || is.na(range) 
      || length(range) != 1 
      || !.is.wholenumber(range) 
      || !(range >= 2 
           && range <= (knnaff$patterns[[knnaff$numPatterns]]$cardinality + 
                  knnaff$patterns[[knnaff$numPatterns - 1]]$cardinality - 1) 
           || range == -1)){
    warning("In treatment number ", i, ": Argument \"range\" not 
            specified correctly. Defaults are applied")
    knnaff$range <- -1
  }else{
    knnaff$range <- range
  }
  if (!is.numeric(k) 
      || is.na(k) 
      || length(k) != 1 
      || !.is.wholenumber(k) 
      || !(k >= 1 
           && k <= (knnaff$patterns[[knnaff$numPatterns]]$cardinality + 
                      knnaff$patterns[[knnaff$numPatterns - 1]]$cardinality) 
           || k == -1)){
    warning("In treatment number ", i, ": Argument \"k\" not specified 
            correctly. Defaults are applied")
    knnaff$k <- -1
  }else{
    knnaff$k <- k
  }
  # Do leave-one-out cross-validation
  knnaff <- knnaff.docv(knnaff)
  
  return (knnaff)
}

knaff.create.structure <- function(data){

  # Elemantary statistics
  dimension <- ncol(data) - 1
  numOfPoints <- nrow(data)
  classNames <- unique(data[,dimension + 1])
  numOfClasses <- length(classNames)
  # Ordering patterns according to their cardinalities
  classCardinalities <- rep(0, numOfClasses)
  for (i in 1:numOfClasses){
    classCardinalities[i] <- nrow(data[data[,dimension + 1] == classNames[i],])
  }
  # Creating pattern templates
  patterns <- as.list("")
  for (i in 1:numOfClasses){
    maxCarIndex <- which.max(classCardinalities)
    # Creating a single template
    pattern.index       <- i
    pattern.points      <- data[data[,dimension + 1] == classNames[maxCarIndex],
                                1:dimension]
    pattern.name        <- classNames[maxCarIndex]
    pattern.cardinality <- classCardinalities[maxCarIndex]
    pattern.votes       <- 0
    pattern <- structure(
      list(index = pattern.index, 
           points = pattern.points, 
           name = pattern.name, 
           cardinality = pattern.cardinality, 
           votes = pattern.votes), 
      .Names = c("index", "points", "name", "cardinality", "votes"))
    # Adding pattern template to the list of patterns
    patterns[[i]] <- pattern
    # Deleting processed pattern
    classCardinalities[maxCarIndex] <- -1    
  }
  # Creating overall structure
  knnaff <- structure(
    list(raw <- data, 
         dimension = dimension, 
         numPatterns = numOfClasses, 
         numPoints = numOfPoints, 
         patterns = patterns, 
         classifiers = list(), 
         numClassifiers = 0, 
         methodAggregation = "majority", 
         range = -1, 
         k = -1), 
    .Names = c("raw", "dimension", "numPatterns", "numPoints", "patterns", 
               "classifiers", "numClassifiers", "methodAggregation", "range", 
               "k"))
  
  return (knnaff)
}

knnaff.docv <- function(knnaff){

  counter <- 1
  # Determining multi-class behaviour
  if (knnaff$methodAggregation == "majority"){
    for (i in 1:(knnaff$numPatterns - 1)){
      for (j in (i + 1):knnaff$numPatterns){
        # Creating a classifier
        classifier.index          <- counter
        classifier.index0         <- i
        classifier.index1         <- j
        classifier.points         <- as.double(t(rbind(knnaff$patterns[[i]]$points, knnaff$patterns[[j]]$points)))
        classifier.cardinalities  <- as.integer(c(knnaff$patterns[[i]]$cardinality, knnaff$patterns[[j]]$cardinality))
        if (knnaff$k < 1 || knnaff$k > (knnaff$patterns[[i]]$cardinality + knnaff$patterns[[j]]$cardinality - 1))
        {
          if (knnaff$range < 2 || knnaff$range > (knnaff$patterns[[i]]$cardinality + knnaff$patterns[[j]]$cardinality - 1)){
            maxk <- 10*( (knnaff$numPoints)^(1/knnaff$dimension) ) + 1
          }else{
            maxk <- knnaff$range
          }
          maxk <- min(maxk, knnaff$patterns[[i]]$cardinality + knnaff$patterns[[j]]$cardinality - 1)
          maxk <- max(maxk, 2)
          classifier.range <- maxk
          classifier.k <- as.integer(.C("KnnAffInvLearnJK", 
                                        classifier.points, 
                                        as.integer(knnaff$dimension), 
                                        classifier.cardinalities, 
                                        as.integer(maxk), 
                                        k=integer(1))$k)
        }else{
          classifier.range <- knnaff$range
          classifier.k <- as.integer(knnaff$k)
        }
        # Adding the classifier to the list of classifiers
        knnaff$classifiers[[counter]] <- structure(
          list(index = classifier.index, 
               index0 = classifier.index0, 
               index1 = classifier.index1, 
               points = classifier.points, 
               cardinalities = classifier.cardinalities, 
               k = classifier.k, 
               range = classifier.range), 
          .Names = c("index", "index0", "index1", "points", "cardinalities", "k", "range"))
        counter <- counter + 1
      }
    }
  }
  if (knnaff$methodAggregation == "sequent"){
    for (i in 1:knnaff$numPatterns){
      anotherClass <- NULL
      for (j in 1:knnaff$numPatterns){
        if (j != i){
          anotherClass <- rbind(anotherClass, knnaff$patterns[[j]]$points)
        }
      }
      classifier.index          <- counter
      classifier.index0         <- i
      classifier.index1         <- -1
      classifier.points         <- as.double(t(rbind(knnaff$patterns[[i]]$points, anotherClass)))
      classifier.cardinalities  <- as.integer(c(knnaff$patterns[[i]]$cardinality, nrow(anotherClass)))
      if (knnaff$k < 1 || knnaff$k > knnaff$numPoints)
      {
        if (knnaff$range < 2 || knnaff$range > (knnaff$numPoints - 1)){
          maxk <- 10*( (knnaff$numPoints)^(1/knnaff$dimension) ) + 1
        }else{
          maxk <- knnaff$range
        }
        maxk <- min(maxk, knnaff$numPoints - 1)
        maxk <- max(maxk, 2)
        classifier.range <- maxk
        classifier.k <- as.integer(.C("KnnAffInvLearnJK", 
                                      classifier.points, 
                                      as.integer(knnaff$dimension), 
                                      classifier.cardinalities, 
                                      as.integer(maxk), 
                                      k=integer(1))$k)
      }else{
        classifier.range <- knnaff$range
        classifier.k <- as.integer(knnaff$k)
      }
      # Adding the classifier to the list of classifiers
      knnaff$classifiers[[counter]] <- structure(
        list(index = classifier.index, 
             index0 = classifier.index0, 
             index1 = classifier.index1, 
             points = classifier.points, 
             cardinalities = classifier.cardinalities, 
             k = classifier.k, 
             range = classifier.range), 
        .Names = c("index", "index0", "index1", "points", "cardinalities", "k", "range"))
      counter <- counter + 1
    }
  }
  
  return (knnaff)
}

knnaff.classify <- function(objects, knnaff){

  # Correct input data
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  # Initialization of the vote array
  votes <- matrix(rep(0, nrow(objects)*knnaff$numPatterns), nrow=nrow(objects), ncol=knnaff$numPatterns)
  for (i in 1:length(knnaff$classifiers)){
    res <- .C("KnnAffInvClassify", 
              as.double(t(objects)), 
              as.integer(nrow(objects)), 
              knnaff$classifiers[[i]]$points, 
              as.integer(knnaff$dimension), 
              knnaff$classifiers[[i]]$cardinalities, 
              knnaff$classifiers[[i]]$k, 
              output=integer(nrow(objects)))$output
    for (j in 1:nrow(objects)){
      if (res[j] == 0){
        votes[j,knnaff$classifiers[[i]]$index0] <- votes[j,knnaff$classifiers[[i]]$index0] + 1
      }else{
        votes[j,knnaff$classifiers[[i]]$index1] <- votes[j,knnaff$classifiers[[i]]$index1] + 1
      }
    }
  }
  # Collect results
  results <- list()
  for (i in 1:nrow(objects)){
    results[[i]] <- knnaff$patterns[[which.max(votes[i,])]]$name
  }
  
  return (results)
}
