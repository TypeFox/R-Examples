ubNCL <-
function(X, Y, k = 3, verbose = TRUE) {
  
  stopifnot(k > 0, class(verbose) == "logical", all(unique(Y) %in% c(0, 1)))
  
  #only numeric features are allowed
  if(any(sapply(X,is.numeric)==FALSE))
    stop("only numeric features are allowed to compute nearest neighbors")
  
  N.orig <- length(Y)
  # If an instance belongs to the majority class and the classification given by its
  # three nearest neighbors contradicts the original class, then it is removed
  data <- ubENN(X, Y, k, verbose)
  X <- data$X
  Y <- data$Y
  
  N <- length(Y)
  i.1 <- which(Y == 1)
  i.0 <- which(Y == 0)
  
  if (length(i.0) == 0) {
    # if there are no 0 obs then don't do anything
    if (verbose) 
      cat("Warning: No remaining majority instances after ENN \n")
    return(list(X = X, Y = Y))
  }
  
  # If an instance belongs to the minority class and its nearest neighbors
  # misclasify it, then the nearest neighbors that belong to the majority class are removed.
  timeRemove <- system.time({
    out.hat <- FNN::knn(train = X, test = X[i.1, ], cl = Y, k = k + 1, prob = TRUE)
    proba.hat <- attr(out.hat, "prob")
    levels(out.hat) <- c(0, 1)
    prob.th <- k/(k+1)
    id.miss <- which((Y[i.1] != out.hat) & (proba.hat >= prob.th))
    if (length(id.miss) == 0) {
      Id <- 1:N
      id2remove <- NULL
    } else {
      out.nn <- attr(out.hat, "nn.index")[, -1]  #the fist column is the point itself
      if (!is.vector(out.nn)) 
        nn.miss <- out.nn[id.miss, ] else nn.miss <- out.nn[id.miss]
      nn.2check <- unique(as.vector(nn.miss))
      id2remove <- nn.2check[which(Y[nn.2check] == 0)]
      Id <- setdiff(1:N, id2remove)
    }
  })
  if (verbose) 
    cat("Number of instances removed from majority class after ENN:", length(id2remove), 
        "\t Time needed:", round(timeRemove[3], digits = 2), "\n")
  
  Id <- sort(Id)
  if (is.vector(X) != TRUE) 
    X = X[Id, ] else X = X[Id]
  Y = Y[Id]
  N.left <- length(Y)
  if (verbose) 
    cat("Number of instances removed from majority class:", N.orig - N.left, "\n")
  
  return(list(X = X, Y = Y))
}
