compclassf.train <- function(dataf, labels, 
                            to.equalize = TRUE, 
                            to.reduce = FALSE, 
                            classifier.type = c("ddalpha", "maxdepth", "knnaff", "lda", "qda"), 
                            ...){
  # Trains the functional componentwise classifier
  # Args:
  #   dataf:  list containing lists (functions) of two vectors of equal length, 
  #           named "args" and "vals": arguments sorted in ascending order and 
  #           corresponding them values respectively
  #   labels: output labels of the functinal observations
  #   other arguments: TODO
  # Returns:
  #   Functional componentwise clasifier
  
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
  
  # Bring to finite dimension
  
  # Pointize
  points <- GetPointsDHB12(dataf, labels, to.equalize, to.reduce)
  # CV
  arg.indices <- getBestSpaceDHB12(points$data, classifier.type, num.chunks=10, ...)
  data <- points$data[,c(arg.indices,ncol(points$data))]
  # Apply chosen classifier to train the data
  if (classifier.type == "ddalpha"){
    classifier <- ddalpha.train(data, separator = "alpha", ...)
  }
  if (classifier.type == "maxdepth"){
    classifier <- ddalpha.train(data, separator = "maxD", ...)
  }
  if (classifier.type == "knnaff"){
    classifier <- knnaff.train(data, i = 0, ...)
  }
  if (classifier.type == "lda"){
    classifier <- lda.train(data, ...)
  }
  if (classifier.type == "qda"){
    classifier <- qda.train(data, ...)
  }
  # Create the eventual output structure
  compclassf <- structure(
    list(dataf = points$dataf, 
      labels = points$labels, 
      adc.method = "equalCover", 
      adc.args = list(instance = "val", numFcn = ncol(points$data) - 1, numDer = 0), 
      adc.transmat = points$transmat, 
      the.args = arg.indices, 
      data = points$data, 
      classifier.type = classifier.type, 
      classifier = classifier), 
    .Names = c("dataf", "labels", "adc.method", "adc.args", "adc.transmat", 
               "the.args",  "data", "classifier.type", "classifier"))
  class(compclassf) <- "compclassf"
  
  return (compclassf)
}

compclassf.classify <- function(objectsf, compclassf, ...){
  # Classifies functions
  # Args:
  #   objectsf: sample to classify, a list containing lists (functions) of 
  #             two vectors of equal length, named "args" and "vals": 
  #             arguments sorted in ascending order and corresponding them 
  #             values respectively
  #   compclassf: functional DDalpha-classifier
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
  if (compclassf$adc.method == "equalCover"){
    if (compclassf$adc.args$instance == "val"){
      input <- getValGrid(objectsf.equalized, 
                          compclassf$adc.args$numFcn, compclassf$adc.args$numDer)
    }
    if (compclassf$adc.args$instance == "avr"){
      input <- getAvrGrid(objectsf.equalized, 
                          compclassf$adc.args$numFcn, compclassf$adc.args$numDer)
    }
    if (!is.null(compclassf$adc.transmat)){
      input <- input%*%compclassf$adc.transmat
    }
  }
  input <- input[,compclassf$the.args]
  # Classify and assign class labels
  if (compclassf$classifier.type == "ddalpha" || compclassf$classifier.type == "maxdepth"){
    output <- ddalpha.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "knnaff"){
    output <- knnaff.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "lda"){
    output <- lda.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "qda"){
    output <- qda.classify(objects = input, compclassf$classifier, ...)
  }
  classes <- list()
  for (i in 1:length(output)){
#    if (is.numeric(output[[i]])){
      classes[[i]] <- compclassf$labels[[ output[[i]] ]]
#    }else{
#      classes[[i]] <- output[[i]]
#    }
  }

  return (classes)
}

print.compclassf <- function(x, ...){
  cat("compclassf:\n")
  cat("\t num.functions = ", length(x$dataf), 
      ", num.patterns = ", length(unique(x$labels)), "\n", sep="")
  #  cat("\t adc.method", x$adc.method, "\"\n", sep="")
  cat("\t adc:", x$adc.args$instance, "; numFcn:", x$adc.args$numFcn, "; numDer:", x$adc.args$numDer, "\"\n", sep="")
  cat("\t adc.transmat", x$adc.transmat, "\"\n", sep="")
  cat("\t classifier.type", x$classifier.type, "\"\n", sep="")
  cat("\t classifier:\n") 
  print(x$classifier)
}

################################################################################
# Functions below are used for intermediate computations                       #
################################################################################

GetPointsDHB12 <- function(dataf, labels, to.equalize=T, to.reduce=F){
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
  # Prepare data
  if (to.equalize){
    num.times = length(dataf[[1]]$args)
    dataf.equalized <- equalize(dataf)
    adc.args = list(instance = "val", 
                    numFcn = num.times, 
                    numDer = 0)
    input <- getValGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
  }else{
    input <- NULL
    for (i in 1:length(dataf)){
      input <- rbind(input, dataf[[i]]$vals)
    }
  }
  transmat <- NULL
  if (to.reduce){# Reduce dimension if needed
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
  }
  # Combine data
  data <- cbind(input, output, deparse.level=0)
  return (list(data = data, dataf = dataf.equalized, labels = names, transmat = transmat))
}

getBestSpaceDHB12 <- function(data, 
                              classifier.type = "ddalpha", 
                              num.chunks = 10, 
                              ...){
  indices.num <- ncol(data) - 1
  indices.avlbl <- rep(TRUE, indices.num)
  indices.best <- c()
  error.last <- nrow(data) + 1
  r <- 0
  while (sum(indices.avlbl) > 0){
    # If this is the first iteration search through all possible pairs
    if (r == 0){
      # Generate all combinations with smallest distance 2
      combinations <- combn((1:indices.num)[indices.avlbl], 2)
      tmp.cmb <- rbind(combinations[-1,], rep(-1000000, ncol(combinations)))
      tmp.cmb <- (tmp.cmb - combinations)==T
      combinations <- combinations[,apply(tmp.cmb, 2, sum)==0]
      # Choose the best combination
      errors <- c()
      for (i in 1:ncol(combinations)){
        cat("r = ", r, ": ", i, "/", ncol(combinations), ".\n", sep="")
        errors <- c(errors, 0)
        # Actually CV
        num.points <- nrow(data)
        indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
        for (j in 1:num.chunks){
          # Determine points to be taken off
          take.off <- (indices.off + j)[(indices.off + j) <= num.points]
          # Apply chosen classifier
          if (classifier.type == "ddalpha"){
            classifier <- ddalpha.train(data[-take.off,c(combinations[,i], indices.num + 1)], separator = "alpha", ...)
            results <- ddalpha.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "maxdepth"){
            classifier <- ddalpha.train(data[-take.off,c(combinations[,i], indices.num + 1)], separator = "maxD", ...)
            results <- ddalpha.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "knnaff"){
            classifier <- knnaff.train(data[-take.off,c(combinations[,i], indices.num + 1)], i = i, ...)
            results <- knnaff.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "lda"){
            classifier <- lda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- lda.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "qda"){
            classifier <- qda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- qda.classify(data[take.off,combinations[,i]], classifier)
          }
          # Collect errors
          errors[i] <- errors[i] + sum(unlist(results) != data[take.off,indices.num + 1])
        }
      }
      # Collect results
      error.last <- min(errors)
      indices.best <- combinations[,which.min(errors)]
      indices.avlbl <- rep(TRUE, indices.num)
      indices.to.dsbl <- unique(c(indices.best, indices.best - 1, indices.best + 1))
      indices.to.dsbl <- indices.to.dsbl[indices.to.dsbl >= 1 && indices.to.dsbl <= indices.num]
      indices.avlbl[indices.to.dsbl] <- FALSE
      r <- 2
      next
    }
    # First, sequential approach
    errors <- c()
    variants <- c()
    for (i in 1:indices.num){
      if (indices.avlbl[i]){
        errors <- c(errors, 0)
        variants <- c(variants, i)
        indices.cur <- c(indices.best, i)
        # Actually CV
        num.points <- nrow(data)
        indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
        for (j in 1:num.chunks){
          # Determine points to be taken off
          take.off <- (indices.off + j)[(indices.off + j) <= num.points]
          # Apply chosen classifier
          if (classifier.type == "ddalpha"){
            classifier <- ddalpha.train(data[-take.off,c(indices.cur, indices.num + 1)], separator = "alpha", ...)
            results <- ddalpha.classify(data[take.off,indices.cur], classifier)
          }
          if (classifier.type == "maxdepth"){
            classifier <- ddalpha.train(data[-take.off,c(indices.cur, indices.num + 1)], separator = "maxD", ...)
            results <- ddalpha.classify(data[take.off,indices.cur], classifier)
          }
          if (classifier.type == "knnaff"){
            classifier <- knnaff.train(data[-take.off,c(indices.cur, indices.num + 1)], i = i, ...)
            results <- knnaff.classify(data[take.off,indices.cur], classifier)
          }
          if (classifier.type == "lda"){
            classifier <- lda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- lda.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "qda"){
            classifier <- qda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- qda.classify(data[take.off,combinations[,i]], classifier)
          }
          # Collect errors
          errors[i] <- errors[i] + sum(unlist(results) != data[take.off,indices.num + 1])
        }
      }
    }
    error.best <- min(errors)
    best.i <- variants[which.min(errors)[1]]
    indices.new <- c(indices.best, best.i)
    # Refinements for r=2, 3 and 4
    if (r %in% 2:3){
      # Define the grid
      if (r == 2){step <- 10}
      if (r == 3){step <- 5}
      grid.one <- c(-(1:step*2), 0, 1:step*2)
      grid <- c()
      for (i in 1:length(indices.new)){
        grid <- c(grid, indices.new[i] + grid.one)
      }
      grid <- unique(grid)
      grid <- sort(grid[(grid >= 1) & (grid <= indices.num)])
      # Generate all combinations with smallest distance 2
      combinations <- combn(grid, r + 1)
      tmp.cmb <- rbind(combinations[-1,], rep(-1000000, ncol(combinations)))
      tmp.cmb <- (tmp.cmb - combinations)==T
      combinations <- combinations[,apply(tmp.cmb, 2, sum)==0]
      # Choose the best combination
      #indices.grid <- (1:indices.num)[indices.avlbl & ((1:induces.num) %in% grid)]
      # Go through the combinations
      errors <- c()
      #combinations <- combn(indices.grid, r + 1)
      for (i in 1:ncol(combinations)){
        cat("r = ", r, ": ", i, "/", ncol(combinations), ".\n", sep="")
        errors <- c(errors, 0)
        # Actually CV
        num.points <- nrow(data)
        indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
        for (j in 1:num.chunks){
          # Determine points to be taken off
          take.off <- (indices.off + j)[(indices.off + j) <= num.points]
          # Apply chosen classifier
          if (classifier.type == "ddalpha"){
            classifier <- ddalpha.train(data[-take.off,c(combinations[,i], indices.num + 1)], separator = "alpha", ...)
            results <- ddalpha.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "maxdepth"){
            classifier <- ddalpha.train(data[-take.off,c(combinations[,i], indices.num + 1)], separator = "maxD", ...)
            results <- ddalpha.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "knnaff"){
            classifier <- knnaff.train(data[-take.off,c(combinations[,i], indices.num + 1)], i = i, ...)
            results <- knnaff.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "lda"){
            classifier <- lda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- lda.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "qda"){
            classifier <- qda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- qda.classify(data[take.off,combinations[,i]], classifier)
          }
          # Collect errors
          errors[i] <- errors[i] + sum(unlist(results) != data[take.off,indices.num + 1])
        }
      }
      error.best <- min(errors)
      indices.cur <- combinations[,which.min(errors)]
    }else{
      indices.cur <- indices.new
    }
    if (error.best < error.last){
      indices.best <- indices.cur
      error.last <- error.best
      indices.avlbl <- rep(TRUE, indices.num)
      indices.to.dsbl <- unique(c(indices.best, indices.best - 1, indices.best + 1))
      indices.to.dsbl <- indices.to.dsbl[indices.to.dsbl >= 1 && indices.to.dsbl <= indices.num]
      indices.avlbl[indices.to.dsbl] <- FALSE
      r <- r + 1
    }else{
      break
    }
  }
  return (indices.best)
}
