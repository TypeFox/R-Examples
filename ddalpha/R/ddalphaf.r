ddalphaf.train <- function(dataf, labels, 
                           adc.args = list(instance = "avr", 
                                           numFcn = -1, 
                                           numDer = -1), 
                           classifier.type = c("ddalpha", "maxdepth", "knnaff", "lda", "qda"), 
                           cv.complete = FALSE, 
                           seed = 0,
                           ...){
  # Trains the functional DDalpha-classifier
  # Args:
  #   dataf:  list containing lists (functions) of two vectors of equal length, 
  #           named "args" and "vals": arguments sorted in ascending order and 
  #           corresponding them values respectively
  #   labels: output labels of the functinal observations
  #   other arguments: TODO
  # Returns:
  #   Functional DDalpha-classifier
  
  # Check "dataf"
  if (!is.list(dataf))
    stop("Argument 'dataf' must be a list")
  for (df in dataf)
    if (!(is.list(df) && length(df) == 2 &&
            !is.null(df$args) && !is.null(df$vals) &&
            is.vector(df$args) && is.vector(df$vals) &&
            is.numeric(df$args) && is.numeric(df$vals) &&
            length(df$args) == length(df$vals) &&
            sort(df$args) == df$args))
      stop("Argument 'dataf' must be a list containing lists (functions) of two vectors of equal length, named 'args' and 'vals': arguments sorted in ascending order and corresponding them values respectively")
  
  # Check "labels"
  if (!(length(dataf)==length(labels) && length(unique(labels)>=2)))
    stop("Argument 'labels' has wrong format")
  
  # Check classifier.type
  classifier.type = match.arg(classifier.type)
  
  # Check "adc.method"
  adc.method = 'equalCover'
  
  if (seed != 0) set.seed(seed)
  
  # Check "adc.args"
  if (!(adc.args$instance %in% c("val", "avr") && 
      ((adc.args$numFcn >= 0 &&  adc.args$numDer >= 0 && (adc.args$numFcn + adc.args$numDer >= 2)) ||
       (adc.args$numFcn == -1 &&  adc.args$numDer == -1))))
    stop("Argument 'adc.args' has wrong format")
  
  # CV
  if (adc.args$numFcn == -1 && adc.args$numDer == -1){
    if (cv.complete){
      res <- getBestSpaceCV(dataf, labels, adc.method, adc.args, 
                            classifier.type, num.chunks=10, ...)
    }else{
      res <- getBestSpace(dataf, labels, adc.method, adc.args,
                          classifier.type, num.chunks=10, ...)
    }
    the.args <- res$args
    num.cv <- res$num.cv
  }else{
    the.args <- adc.args
    num.cv <- 0
  }
  # Pointize
  points <- GetPoints(dataf, labels, adc.method, the.args)
  # Apply chosen classifier to train the data
  if (classifier.type == "ddalpha"){
    classifier <- ddalpha.train(points$data, seed = seed, ...)
  }
  if (classifier.type == "maxdepth"){
    classifier <- ddalpha.train(points$data, separator = "maxD", seed = seed, ...)
  }
  if (classifier.type == "knnaff"){
    classifier <- knnaff.train(points$data, i = 0, ...)
  }
  if (classifier.type == "lda"){
    classifier <- lda.train(points$data, ...)
  }
  if (classifier.type == "qda"){
    classifier <- qda.train(points$data, ...)
  }
  # Create the eventual output structure
  ddalphaf <- structure(
    list(dataf = points$dataf, 
      labels = points$labels, 
      adc.method = adc.method, 
      adc.args = the.args, 
      adc.num.cv = num.cv, 
      adc.transmat = points$adc.transmat, 
      data = points$data, 
      classifier.type = classifier.type, 
      classifier = classifier), 
    .Names = c("dataf", "labels", "adc.method", "adc.args", "adc.num.cv", 
               "adc.transmat", "data", "classifier.type", "classifier"))
  class(ddalphaf) <- "ddalphaf"
  
  return (ddalphaf)
}

ddalphaf.classify <- function(objectsf, ddalphaf, ...){
  # Classifies functions
  # Args:
  #   objectsf: sample to classify, a list containing lists (functions) of 
  #             two vectors of equal length, named "args" and "vals": 
  #             arguments sorted in ascending order and corresponding them 
  #             values respectively
  #   ddalphaf: functional DDalpha-classifier
  # Returns:
  #   List of labels assigned to the functions from "objectsf"

  # Check "objectsf"
  if (!is.list(objectsf))
    stop("Argument 'objectsf' must be a list")
  if (!is.null(objectsf$args)){
    objectsf = list(objectsf) # there was a single element 
  }
  for (df in objectsf)
    if (!(is.list(df) && length(df) == 2 &&
            !is.null(df$args) && !is.null(df$vals) &&
            is.vector(df$args) && is.vector(df$vals) &&
            is.numeric(df$args) && is.numeric(df$vals) &&
            length(df$args) == length(df$vals) &&
            sort(df$args) == df$args))
      stop("Argument 'objectsf' must be a list containing lists (functions) of two vectors of equal length, named 'args' and 'vals': arguments sorted in ascending order and corresponding them values respectively")

  # Prepare to multivariate classification
  objectsf.equalized <- equalize(objectsf)
  if (ddalphaf$adc.method == "equalCover"){
    if (ddalphaf$adc.args$instance == "val"){
      input <- getValGrid(objectsf.equalized, 
                          ddalphaf$adc.args$numFcn, ddalphaf$adc.args$numDer)
    }
    if (ddalphaf$adc.args$instance == "avr"){
      input <- getAvrGrid(objectsf.equalized, 
                          ddalphaf$adc.args$numFcn, ddalphaf$adc.args$numDer)
    }
    if (!is.null(ddalphaf$adc.transmat)){
      input <- input%*%ddalphaf$adc.transmat
    }
  }
  # Classify and assign class labels
  if (ddalphaf$classifier.type == "ddalpha" || ddalphaf$classifier.type == "maxdepth"){
    output <- ddalpha.classify(objects = input, ddalphaf$classifier, ...)
  }
  if (ddalphaf$classifier.type == "knnaff"){
    output <- knnaff.classify(objects = input, ddalphaf$classifier, ...)
  }
  if (ddalphaf$classifier.type == "lda"){
    output <- lda.classify(objects = input, ddalphaf$classifier, ...)
  }
  if (ddalphaf$classifier.type == "qda"){
    output <- qda.classify(objects = input, ddalphaf$classifier, ...)
  }
  classes <- list()
  for (i in 1:length(output)){
#    if (is.numeric(output[[i]])){
      classes[[i]] <- ddalphaf$labels[[ output[[i]] ]]
#    }else{
#      classes[[i]] <- output[[i]]
#    }
  }

  return (classes)
}

is.in.convexf <- function(objectsf, dataf, cardinalities, 
                          adc.method = "equalCover", 
                          adc.args = list(instance = "val", 
                                          numFcn = 5, 
                                          numDer = 5), seed = 0){
  # Checks if the function(s) lie(s) inside convex hulls of the 
  # functions from the sample in the projection space
  # Args:
  #   objectsf:      list containing lists (functions) of two vectors of equal 
  #                  length, named "args" and "vals": arguments sorted in 
  #                  ascending order and corresponding them values 
  #                  respectively. These functions are supposed to be checked 
  #                  for 'outsiderness'
  #   dataf:         list containing lists (functions) of two vectors of equal 
  #                  length, named "args" and "vals": arguments sorted in 
  #                  ascending order and corresponding them values respectively
  #   cardinalities: cardinalities of the classes in "dataf"
  #   other arguments: TODO
  # Returns:
  #   Vector with 1s for those lying inside of at least one of the convex hulls 
  #   of the classes and 0s for those lying beyond them
  
  # Data-consistency checks
  # TODO
  
  # Project "objectsf" into a multivariate space
  objectsf.equalized <- equalize(objectsf)
  if (adc.method == "equalCover"){
    if (adc.args$instance == "val"){
      objects <- getValGrid(objectsf.equalized, adc.args$numFcn, adc.args$numDer)
    }
    if (adc.args$instance == "avr"){
      objects <- getAvrGrid(objectsf.equalized, adc.args$numFcn, adc.args$numDer)
    }
  }
  # Project "dataf" into a multivariate space
  dataf.equalized <- equalize(dataf)
  if (adc.method == "equalCover"){
    if (adc.args$instance == "val"){
      data <- getValGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
    }
    if (adc.args$instance == "avr"){
      data <- getAvrGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
    }    
  }
  in.convex <- is.in.convex(objects, data, cardinalities, seed)
  
  return (in.convex)
}

print.ddalphaf <- function(x, ...){
  cat("ddalphaf:\n")
  cat("\t num.functions = ", length(x$dataf), 
  ", num.patterns = ", length(unique(x$labels)), "\n", sep="")
#  cat("\t adc.method", x$adc.method, "\"\n", sep="")
  cat("\t adc:", x$adc.args$instance, "; numFcn:", x$adc.args$numFcn, "; numDer:", x$adc.args$numDer, "\"\n", sep="")
  cat("\t adc.num.cv", x$adc.num.cv, "\"\n", sep="")
  cat("\t adc.transmat", x$adc.transmat, "\"\n", sep="")
  cat("\t classifier.type", x$classifier.type, "\"\n", sep="")
  cat("\t classifier:\n") 
  print(x$classifier)
}

################################################################################
# Functions below are used for intermediate computations                       #
################################################################################

equalize <- function(dataf){
  # 1. Adjusts the data to have equal (the largest) argument interval
  # 2. Calclates - numerically - derivative
  # Args:
  #   dataf: list containing lists (functions) of two vectors 
  #          of equal length, named "args" and "vals": arguments 
  #          sorted in ascending order and corresponding them 
  #          values respectively
  # Returns:
  #   The list of lists of the same structure, 'equalized', 
  #   contating derivatives as 3rd vector named "der1"
  
  # Check whether every function contains fields "args" and "vals", 
  # whether they are numerical and of equal length, have no NAs or 
  # ties and are sorted in ascending order
  # TODO
  
  # 1.
  # Get argument bounds
  min <- Inf
  max <- -Inf
  for (i in 1:length(dataf)){
    if (dataf[[i]]$args[1] < min){min <- dataf[[i]]$args[1]}
    if (dataf[[i]]$args[length(dataf[[i]]$args)] > max){
      max <- dataf[[i]]$args[length(dataf[[i]]$args)]
    }
  }
  # and apply them to equalize functions timely
  for (i in 1:length(dataf)){
    if (dataf[[i]]$args[1] > min){
      dataf[[i]]$args <- c(min, dataf[[i]]$args)
      dataf[[i]]$vals <- c(dataf[[i]]$vals[1], dataf[[i]]$vals)
    }
    if (dataf[[i]]$args[length(dataf[[i]]$args)] < max){
      dataf[[i]]$args <- c(dataf[[i]]$args, max)
      dataf[[i]]$vals <- c(dataf[[i]]$vals, 
                           dataf[[i]]$vals[length(dataf[[i]]$vals)])
    }
    # Computational trick - add "-1" to the "left"
    #dataf[[i]]$args <- c(min - 1, dataf[[i]]$args)
    #dataf[[i]]$vals <- c(dataf[[i]]$vals[1], dataf[[i]]$vals)
  }
  
  # 2.
  for (i in 1:length(dataf)){
    dataf[[i]] <- derive(dataf[[i]])
  }
  
  return (dataf)
}

derive <- function(fcn){
  # Adds 1st derivative to the function: a vector named "der1"
  # Args:
  #   fcn: function, a list of two vectors of equal length, named "args" 
  #   and "vals": arguments sorted in ascending order and corresponding 
  #   them values respectively
  # Returns:
  #   The list of of the same structure contating derivative as 3rd 
  #   vector named "der1"
  
  fcn$der1 <- rep(0, length(fcn$args))
  fcn$der1[1] = 0
  for (i in 2:length(fcn$der1)){
    fcn$der1[i] <- (fcn$vals[i] - fcn$vals[i - 1])/
      (fcn$args[i] - fcn$args[i - 1])
  }
  
  return (fcn)
}

getValue <- function(fcn, arg, fcnInstance){
  # Gets the value of the function or its derivative for the given argument 
  # value
  # Args:
  #   fcn:         function, a list of vectors of equal length, named "args" 
  #                (arguments), "vals" (function values) [and it's 
  #                derivatives of order "I", named derI]; arguments (and 
  #                corresponding values) sorted in ascending order
  #   arg:         argument value at which the function (derivative) value 
  #                is to be taken
  #   fcnInstance: inctance to be evaluated; "vals" for the function values, 
  #                "der1" for the first derivative
  # Returns:
  #   Value of the function (derivative)
  
  # Check "arg", evt. "fcnInstance"
  # TODO
  
  # Find corresponding interval
  index <- 2
  while (arg > fcn$args[index]){
    index <- index + 1
  }
  # Get the function value(s) by linear approximation
  if (fcnInstance == "vals"){
    value <- fcn$vals[index - 1] + 
      (fcn$vals[index] - fcn$vals[index - 1])*
      ((arg - fcn$args[index - 1])/(fcn$args[index] - fcn$args[index - 1]))
  }
  if (fcnInstance == "der1"){
    value <- fcn$der1[index]
  }
  
  return (value)
}

getAvrValue <- function(fcn, argFrom, argTo, fcnInstance){
  # Gets the average value of the function or its derivative on the given 
  # interval
  # Args:
  #   fcn:         function, a list of vectors of equal length, named "args" 
  #                (arguments), "vals" (function values) [and it's 
  #                derivatives of order "I", named derI]; arguments (and 
  #                corresponding values) sorted in ascending order
  #   argFrom:     argument value from which the function (derivative) value 
  #                is to be averaged
  #   argTo:       argument value to which the function (derivative) value 
  #                is to be averaged
  #   fcnInstance: inctance to be evaluated; "vals" for the function values, 
  #                "der1" for the first derivative
  # Returns:
  #   Average value of the function (derivative) on the interval 
  #   (argFrom, argTo)
  
  # Check "argFrom" and "argTo", evt. "fcnInstance"
  # TODO
  
  # Define 'from' and 'to' interval
  indexFrom <- 2
  while (argFrom > fcn$args[indexFrom]){
    indexFrom <- indexFrom + 1
  }
  indexTo <- 2
  while (argTo > fcn$args[indexTo]){
    indexTo <- indexTo + 1
  }
  average <- 0
  valTo <- getValue(fcn, argTo, fcnInstance)
  # Integrate
  curArgFrom <- argFrom
  if (fcnInstance == "vals"){
    valFrom <- getValue(fcn, curArgFrom, "vals")
    while(indexFrom < indexTo){
      average <- average + (valFrom + fcn$vals[indexFrom])*
        (fcn$args[indexFrom] - curArgFrom)/2
      valFrom <- fcn$vals[indexFrom]
      curArgFrom <- fcn$args[indexFrom]
      indexFrom <- indexFrom + 1
    }
    average <- average + (valFrom + valTo)*(argTo - curArgFrom)/2
  }
  if (fcnInstance == "der1"){
    while(indexFrom < indexTo){
      average <- average + (fcn$der1[indexFrom])*
        (fcn$args[indexFrom] - curArgFrom)
      curArgFrom <- fcn$args[indexFrom]
      indexFrom <- indexFrom + 1
    }
    average <- average + valTo*(argTo - curArgFrom)
  }
  average <- average/(argTo - argFrom)
  
  return (average)
}

getValGrid <- function(dataf, numFcn, numDer){
  # Represents a function sample as a multidimensional (d="numFcn"+"numDer") 
  # one averaging for that each function and it derivative on "numFcn" 
  # (resp. "numDer") equal nonoverlapping covering intervals
  # Args:
  #   dataf:  list containing lists (functions) of vectors of equal length, 
  #           first two named "args" and "vals" are arguments sorted in 
  #           ascending order and having same bounds for all functions and 
  #           corresponding them values respectively
  #   numFcn: number of function intervals
  #   numDer: number of first-derivative intervals
  # Returns:
  #   Matrix - a multidimensional presentation of the functional sample
  
  # Get argument bounds ("dataf" is equalized)
  min <- dataf[[1]]$args[1]
  max <- dataf[[1]]$args[length(dataf[[1]]$args)]
  # Get argument grid
  args <- dataf[[1]]$args
  argsFcn <- min + 0:numFcn*(max - min)/(numFcn - 1)
  argsDer <- min + 0:numDer*(max - min)/(numDer - 1)
  # Get function/derivative grid
  fcnGrid <- matrix(nrow = length(dataf), ncol = numFcn)
  derGrid <- matrix(nrow = length(dataf), ncol = numDer)
  if (numFcn > 0){
    for (i in 1:length(dataf)){
      # Set merging algorithm (Novikoff, Piter)
      cArgs <- 1
      cArgsFcn <- 1
      fcnGrid[i,1] <- dataf[[i]]$vals[1]
      while (cArgsFcn != numFcn){
#        print(argsFcn)
#        print(fcnGrid[i,])
#        cat(cArgs, " and ", cArgsFcn, "\n")
#        cat(args[cArgs + 1], " and ", argsFcn[cArgsFcn + 1], "\n")
        if (args[cArgs + 1] < argsFcn[cArgsFcn + 1]){
          cArgs <- cArgs + 1
        }else{
          nextArg <- argsFcn[cArgsFcn + 1]
          fcnGrid[i,cArgsFcn + 1] <- dataf[[i]]$vals[cArgs] + (nextArg - args[cArgs])*dataf[[i]]$der1[cArgs + 1]
          if (args[cArgs + 1] == argsFcn[cArgsFcn + 1]){
            cArgs <- cArgs + 1
          }
          cArgsFcn <- cArgsFcn + 1
        }
      }
    }
  }
  if (numDer > 0){
    for (i in 1:length(dataf)){
      # Again, set merging algorithm (Novikoff, Piter)
      cArgs <- 1
      cArgsDer <- 1
      derGrid[1] <- dataf[[i]]$ders[2]
      while (cArgsDer != numDer){
#        print(argsDer)
#        print(derGrid[i,])
#        cat(cArgs, " and ", cArgsDer, "\n")
#        cat(args[cArgs + 1], " and ", argsDer[cArgsDer + 1], "\n")
        if (args[cArgs + 1] < argsDer[cArgsDer + 1]){
          cArgs <- cArgs + 1
        }else{
          derGrid[i,cArgsDer + 1] <- dataf[[i]]$der1[cArgs + 1]
          if (args[cArgs + 1] == argsDer[cArgsDer + 1]){
            cArgs <- cArgs + 1
          }
          cArgsDer <- cArgsDer + 1
        }
      }
    }
  }
  mvX <- cbind(fcnGrid, derGrid)
  return (mvX)
}

getAvrGrid <- function(dataf, numFcn, numDer){
  # Represents a function sample as a multidimensional (d="numFcn"+"numDer") 
  # one averaging for that each function and it derivative on "numFcn" 
  # (resp. "numDer") equal nonoverlapping covering intervals
  # Args:
  #   dataf:  list containing lists (functions) of vectors of equal length, 
  #           first two named "args" and "vals" are arguments sorted in 
  #           ascending order and having same bounds for all functions and 
  #           corresponding them values respectively
  #   numFcn: number of function intervals
  #   numDer: number of first-derivative intervals
  # Returns:
  #   Matrix - a multidimensional presentation of the functional sample
  
  # Get argument bounds ("dataf" is equalized)
  min <- dataf[[1]]$args[1]
  max <- dataf[[1]]$args[length(dataf[[1]]$args)]
  # Get argument grid
  args <- dataf[[1]]$args
  argsFcn <- min + 0:numFcn*(max - min)/numFcn
  argsDer <- min + 0:numDer*(max - min)/numDer
  # Get function/derivative grid
  fcnGrid <- matrix(nrow = length(dataf), ncol = numFcn)
  derGrid <- matrix(nrow = length(dataf), ncol = numDer)
  if (numFcn > 0){
    for (i in 1:length(dataf)){
      # Set merging algorithm (Novikoff, Piter)
      cArgs <- 1
      cArgsFcn <- 1
      curArg <- min
      curFcn <- dataf[[i]]$vals[1]
      curAvr <- 0
      while (cArgsFcn != numFcn + 1){
#        print(argsFcn)
#        print(fcnGrid[i,])
#        cat(cArgs, " and ", cArgsFcn, "\n")
#        cat(args[cArgs + 1], " and ", argsFcn[cArgsFcn + 1], "\n")
        if (args[cArgs + 1] < argsFcn[cArgsFcn + 1]){
          nextArg <- args[cArgs + 1]
          nextFcn <- dataf[[i]]$vals[cArgs + 1]
          curAvr <- curAvr + (nextArg - curArg)*(nextFcn + curFcn)/2
          cArgs <- cArgs + 1
        }else{
          nextArg <- argsFcn[cArgsFcn + 1]
          nextFcn <- dataf[[i]]$vals[cArgs] + (nextArg - args[cArgs])*dataf[[i]]$der1[cArgs + 1]
          fcnGrid[i,cArgsFcn] <- curAvr + (nextArg - curArg)*(nextFcn + curFcn)/2
          curAvr <- 0
          if (args[cArgs + 1] == argsFcn[cArgsFcn + 1]){
            cArgs <- cArgs + 1
          }
          cArgsFcn <- cArgsFcn + 1
        }
        curArg <- nextArg
        curFcn <- nextFcn
      }
    }
  }
  fcnGrid <- fcnGrid/(argsFcn[2] - argsFcn[1])
  if (numDer > 0){
    for (i in 1:length(dataf)){
      # Again, set merging algorithm (Novikoff, Piter)
      cArgs <- 1
      cArgsDer <- 1
      curArg <- min
      curAvr <- 0
      while (cArgsDer != numDer + 1){
#        print(argsDer)
#        print(derGrid[i,])
#        cat(cArgs, " and ", cArgsDer, "\n")
#        cat(args[cArgs + 1], " and ", argsDer[cArgsDer + 1], "\n")
        if (args[cArgs + 1] < argsDer[cArgsDer + 1]){
          nextArg <- args[cArgs + 1]
          curAvr <- curAvr + (nextArg - curArg)*dataf[[i]]$der1[cArgs + 1]
          cArgs <- cArgs + 1
        }else{
          nextArg <- argsDer[cArgsDer + 1]
          derGrid[i,cArgsDer] <- curAvr + (nextArg - curArg)*dataf[[i]]$der1[cArgs + 1]
          curAvr <- 0
          if (args[cArgs + 1] == argsDer[cArgsDer + 1]){
            cArgs <- cArgs + 1
          }
          cArgsDer <- cArgsDer + 1
        }
        curArg <- nextArg
      }
    }
  }
  derGrid <- derGrid/(argsDer[2] - argsDer[1])
  mvX <- cbind(fcnGrid, derGrid)
  return (mvX)
}

getVapnikBound <- function(points, dim = NULL){
  n <- nrow(points)
  d <- ncol(points) - 1
  lda <- lda.train(points)
  result <- lda.classify(points[,1:d], lda)
  empRisk <- sum(result != points[,d + 1])/n
  # Calculate the deviation from the empirical risk
  nu <- 1/n
  C <- 0
  for (k in 0:d){
    C <- C + 2*choose(n - 1, k)
  }
  epsilon <- sqrt( (log(C) - log(nu)) / (2*n) )
  
  return (empRisk + epsilon)
}

GetPoints <- function(dataf, labels, adc.method, adc.args){
  # Numerize labels
  names <- unique(labels)
  output <- rep(0, length(labels))
  for (i in 1:length(labels)){
    for (j in 1:length(names)){
      if (labels[[i]] == names[[j]]){
        output[i] = j
        break
      }
    }
  }
  # Pointize data
  dataf.equalized <- equalize(dataf)
  if (adc.method == "equalCover"){
    if (adc.args$instance == "val"){
      input <- getValGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
    }
    if (adc.args$instance == "avr"){
      input <- getAvrGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
    }
    # Reduce dimension if needed
    princomps <- NULL
    newDim <- ncol(input)
    for (i in 1:length(names)){
      classi <- input[output == i,1:ncol(input)]
      princompsi <- prcomp(x=classi, tol=sqrt(.Machine$double.eps))
      #print(princompsi$sdev)
      newDimi <- sum(princompsi$sdev > sqrt(.Machine$double.eps))
      if (newDimi < newDim){
        newDim <- newDimi
        princomps <- princompsi
      }
    }
    #print(newDim)
    transmat <- NULL
    if (newDim < ncol(input)){
      transmat <- matrix(as.vector(princomps$rotation[,1:newDim]), ncol=newDim)
      input <- input%*%transmat
    }
    # Combine data
    data <- cbind(input, output, deparse.level=0)
  }
  
  return (list(data = data, 
               dataf = dataf.equalized, 
               labels = names, 
               adc.transmat = transmat))
}

GetPointsAll <- function(dataf, labels, adc.method = "equalCover", 
                         adc.args = list(instance = "avr", 
                                         numFcn = -1, 
                                         numDer = -1)){
  # Numerize labels
  names <- unique(labels)
  output <- rep(0, length(labels))
  for (i in 1:length(labels)){
    for (j in 1:length(names)){
      if (labels[[i]] == names[[j]]){
        output[i] = j
        break
      }
    }
  }
  # Prepare values
  numMax <- ceiling(length(dataf[[1]]$args)/2 + sqrt(.Machine$double.eps))
  min <- dataf[[1]]$args[1]
  max <- dataf[[1]]$args[length(dataf[[1]]$args)]
  args <- dataf[[1]]$args
  fcnsAll <- list("")
  dersAll <- list("")
  dataf.equalized <- equalize(dataf)
  # Generate all one-type-argument ("fcn" or "der") for all num(Fcn,Der)-values
  for (numFcn in 1:numMax){
    numDer <- numFcn
    argsFcn <- min + 0:numFcn*(max - min)/numFcn
    argsDer <- min + 0:numDer*(max - min)/numDer
    # Get function/derivative grid
    fcnGrid <- matrix(nrow = length(dataf.equalized), ncol = numFcn)
    derGrid <- matrix(nrow = length(dataf.equalized), ncol = numDer)
    for (i in 1:length(dataf.equalized)){
      # Set merging algorithm (Novikoff, Piter)
      cArgs <- 1
      cArgsFcn <- 1
      curArg <- min
      curFcn <- dataf.equalized[[i]]$vals[1]
      curAvr <- 0
      while (cArgsFcn != numFcn + 1){
        if (args[cArgs + 1] < argsFcn[cArgsFcn + 1]){
          nextArg <- args[cArgs + 1]
          nextFcn <- dataf.equalized[[i]]$vals[cArgs + 1]
          curAvr <- curAvr + (nextArg - curArg)*(nextFcn + curFcn)/2
          cArgs <- cArgs + 1
        }else{
          nextArg <- argsFcn[cArgsFcn + 1]
          nextFcn <- dataf.equalized[[i]]$vals[cArgs] + (nextArg - args[cArgs])*dataf.equalized[[i]]$der1[cArgs + 1]
          fcnGrid[i,cArgsFcn] <- curAvr + (nextArg - curArg)*(nextFcn + curFcn)/2
          curAvr <- 0
          if (args[cArgs + 1] == argsFcn[cArgsFcn + 1]){
            cArgs <- cArgs + 1
          }
          cArgsFcn <- cArgsFcn + 1
        }
        curArg <- nextArg
        curFcn <- nextFcn
      }
    }
    fcnsAll[[numFcn]] <- fcnGrid/(argsFcn[2] - argsFcn[1])
    for (i in 1:length(dataf.equalized)){
      # Again, set merging algorithm (Novikoff, Piter)
      cArgs <- 1
      cArgsDer <- 1
      curArg <- min
      curAvr <- 0
      while (cArgsDer != numDer + 1){
        if (args[cArgs + 1] < argsDer[cArgsDer + 1]){
          nextArg <- args[cArgs + 1]
          curAvr <- curAvr + (nextArg - curArg)*dataf.equalized[[i]]$der1[cArgs + 1]
          cArgs <- cArgs + 1
        }else{
          nextArg <- argsDer[cArgsDer + 1]
          derGrid[i,cArgsDer] <- curAvr + (nextArg - curArg)*dataf.equalized[[i]]$der1[cArgs + 1]
          curAvr <- 0
          if (args[cArgs + 1] == argsDer[cArgsDer + 1]){
            cArgs <- cArgs + 1
          }
          cArgsDer <- cArgsDer + 1
        }
        curArg <- nextArg
      }
    }
    dersAll[[numDer]] <- derGrid/(argsDer[2] - argsDer[1])
  }
  pointsAll <- list("")
  counter <- 1
  # Construct the spaces, reducing dimension if needed
  for (dim in 2:numMax){
    for(nDer in 0:dim){
      nFcn <- dim - nDer
      tmp.args <- list(instance=adc.args$instance, numFcn=nFcn, numDer=nDer)
      if (nFcn == 0){
        input <- dersAll[[nDer]]
      }else{
        if (nDer == 0){
          input <- fcnsAll[[nFcn]]
        }else{
          input <- cbind(fcnsAll[[nFcn]], dersAll[[nDer]])
        }
      }
      # Reduce dimension if needed
      princomps <- NULL
      newDim <- ncol(input)
      for (i in 1:length(names)){
        classi <- input[output == i,1:ncol(input)]
        princompsi <- prcomp(x=classi, tol=sqrt(.Machine$double.eps))
        newDimi <- sum(princompsi$sdev > sqrt(.Machine$double.eps))
        if (newDimi < newDim){
          newDim <- newDimi
          princomps <- princompsi
        }
      }
      transmat <- NULL
      if (newDim < ncol(input)){
        transmat <- matrix(as.vector(princomps$rotation[,1:newDim]), ncol=newDim)
        input <- input%*%transmat
      }
      # Combine data
      tmp.points <- cbind(input, output, deparse.level=0)
      # Save to the list
      pointsAll[[counter]] <- list(data = tmp.points, 
                                   adc.args = tmp.args, 
                                   adc.transmat = transmat)
      counter <- counter + 1
    }
  }
  return (pointsAll)
}

getBestSpace <- function(dataf, labels, adc.method = "equalCover", 
                         adc.args = list(instance = "avr", 
                                         numFcn = -1, 
                                         numDer = -1), 
                         classifier.type = "ddalpha", 
                         num.chunks = 10, 
                         ...){
  # First, get Vapnik bounds for all plausible spaces
  numMax <- ceiling(length(dataf[[1]]$args)/2 + sqrt(.Machine$double.eps))
  numTries <- numMax * (numMax + 1) / 2 + numMax + 1 - 3
  Vapnik.bounds <- rep(Inf, numTries)
  numsFcn <- rep(Inf, numTries)
  numsDer <- rep(Inf, numTries)
  curTry <- 1
  pointsAll <- GetPointsAll(dataf, labels, adc.method, adc.args)  
  for (i in 1:length(pointsAll)){
    tmp.args <- pointsAll[[i]]$adc.args
    tmp.points <- pointsAll[[i]]$data
    Vapnik.bounds[curTry] <- getVapnikBound(tmp.points, ncol(tmp.points) - 1)
    numsFcn[curTry] <- pointsAll[[i]]$adc.args$numFcn
    numsDer[curTry] <- pointsAll[[i]]$adc.args$numDer
    curTry <- curTry + 1
  }
  # Second, get 5 best (i.e. lowest) Vapnik bounds
  #   etalons <- sort(Vapnik.bounds)[1:5]
  #   best.indices <- c()
  #   while (length(best.indices) < 5 && length(etalons) > 0){
  #     best.indices <- c(best.indices, which(Vapnik.bounds == etalons[1]))
  #     etalons <- etalons[-1]
  #   }
  #   best.indices <- best.indices[1:5]
  best.indices <- which(Vapnik.bounds %in% sort(Vapnik.bounds)[1:5])
  # Third, cross-validate over these best spaces
  errors <- rep(0, length(best.indices))
  for (i in 1:length(best.indices)){
    tmp.args <- adc.args
    tmp.args$numFcn <- numsFcn[best.indices[i]]
    tmp.args$numDer <- numsDer[best.indices[i]]
    points.all <- pointsAll[[best.indices[i]]]$data
    d <- ncol(points.all) - 1
    # Actually CV
    num.points <- nrow(points.all)
    indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
    for (j in 1:num.chunks){
      # Determine points to be taken off
      take.off <- (indices.off + j)[(indices.off + j) <= num.points]
      # Apply chosen classifier
      if (classifier.type == "ddalpha"){
        classifier <- ddalpha.train(points.all[-take.off,], ...)
        results <- ddalpha.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "maxdepth"){
        classifier <- ddalpha.train(points.all[-take.off,], separator = "maxD", ...)
        results <- ddalpha.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "knnaff"){
        classifier <- knnaff.train(points.all[-take.off,], i = i, ...)
        results <- knnaff.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "lda"){
        classifier <- lda.train(points.all[-take.off,], ...)
        results <- lda.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "qda"){
        classifier <- qda.train(points.all[-take.off,], ...)
        results <- qda.classify(points.all[take.off,1:d], classifier)
      }
      # Collect errors
      errors[i] <- errors[i] + sum(
        unlist(results) != points.all
        [take.off,d + 1])
    }
  }
  best.i <- which.min(errors)
  new.args <- adc.args
  new.args$numFcn <- numsFcn[best.indices[best.i]]
  new.args$numDer <- numsDer[best.indices[best.i]]
  return (list(args = new.args, num.cv = length(best.indices)))
}

getBestSpaceCV <- function(dataf, labels, adc.method = "equalCover", 
                           adc.args = list(instance = "val", 
                                           numFcn = -1, 
                                           numDer = -1), 
                           classifier.type = "ddalpha", 
                           num.chunks = 10, 
                         ...){
  numMax <- ceiling(length(dataf[[1]]$args)/2 + sqrt(.Machine$double.eps))
  numTries <- numMax * (numMax + 1) / 2 + numMax + 1 - 3
  curTry <- 1
  pointsAll <- GetPointsAll(dataf, labels, adc.method, adc.args)
  errors <- rep(0, length(pointsAll))
  for (i in 1:length(pointsAll)){
    points.all <- pointsAll[[i]]$data
    d <- ncol(points.all) - 1
    # Actually CV
    num.points <- nrow(points.all)
    indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
    for (j in 1:num.chunks){
      # Determine points to be taken off
      take.off <- (indices.off + j)[(indices.off + j) <= num.points]
      # Apply chosen classifier
      if (classifier.type == "ddalpha"){
        classifier <- ddalpha.train(points.all[-take.off,], ...)
        results <- ddalpha.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "maxdepth"){
        classifier <- ddalpha.train(points.all[-take.off,], separator = "maxD", ...)
        results <- ddalpha.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "knnaff"){
        classifier <- knnaff.train(points.all[-take.off,], i = i, ...)
        results <- knnaff.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "lda"){
        classifier <- lda.train(points.all[-take.off,], ...)
        results <- lda.classify(points.all[take.off,1:d], classifier)
      }
      if (classifier.type == "qda"){
        classifier <- qda.train(points.all[-take.off,], ...)
        results <- qda.classify(points.all[take.off,1:d], classifier)
      }
      # Collect errors
      errors[i] <- errors[i] + sum(
        unlist(results) != points.all
        [take.off,d + 1])
    }
  }
  best.i <- which.min(errors)
  new.args <- adc.args
  new.args$numFcn <- pointsAll[[best.i]]$adc.args$numFcn
  new.args$numDer <- pointsAll[[best.i]]$adc.args$numDer
  return (list(args = new.args, num.cv = length(pointsAll)))
}
