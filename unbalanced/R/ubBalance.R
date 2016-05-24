ubBalance <-
function(X, Y, type = "ubSMOTE", positive = 1, 
                      percOver = 200, percUnder = 200, k = 5, 
                      perc = 50, method = "percPos", w = NULL, 
                      verbose = FALSE) {
  
  if (any(is.na(Y))) 
    stop("Y has NAs")
  
  if (!is.factor(Y)) 
    stop("Y must be a factor")
  
  lev <- levels(Y)
  if (length(lev) != 2) 
    stop("Y must be a binary factor variable")
  
  # transform the output in the range {0, 1}
  Y <- factor(Y == positive, levels = c(FALSE, TRUE), labels = c(0, 1))
  
  
  if (length(type) > 1) 
    stop("balance type does not support multiple selection")
  
  N.0 <- length(which(Y == 0))
  N.1 <- length(which(Y == 1))
  if (N.0 == 0) {
    cat("Warning: No negative instances, skip balance \n")
    return(list(X = X, Y = Y))
  }
  
  if (N.1 == 0) {
    cat("Warning: No positive instances, skip balance \n")
    return(list(X = X, Y = Y))
  }
  
  if (N.0 == N.1) {
    cat("Warning: equal number of positives and negatives, skip balance \n")
    return(list(X = X, Y = Y))
  }
  
  if (N.0 < N.1) 
    stop(positive, " class is not the minority class")
  
  data <- NULL
  
  if (type == "ubOver") 			data <- ubOver(X, Y, k, verbose)
  if (type == "ubUnder") 	    data <- ubUnder(X, Y, perc, method, w)
  if (type == "ubSMOTE") 	    data <- ubSMOTE(X, Y, percOver, k, percUnder, verbose)
  if (type == "ubOSS") 	      data <- ubOSS(X, Y, verbose)
  if (type == "ubCNN") 	      data <- ubCNN(X, Y, k, verbose)
  if (type == "ubENN") 	      data <- ubENN(X, Y, k, verbose)
  if (type == "ubNCL") 	      data <- ubNCL(X, Y, k, verbose)
  if (type == "ubTomek") 	    data <- ubTomek(X, Y, verbose)
  
  if (is.null(data)) 
    stop("technique", type, " not supported")
  
  X <- data$X
  Y <- data$Y
  id.rm <- data$id.rm
  if (is.null(id.rm)) 
    id.rm <- NA
  
  N <- length(Y)
  
  # Id <- sample(1:N)
  # if (!is.vector(X)) 
  # X = X[Id, ] 
  # else {
  # # is.vector
  # X = X[Id]
  # if (any(is.na(X))) 
  # cat("WARNINGS: vector has NAs \n")
  # if (all(X == X[1])) 
  # cat("WARNINGS: constant vector after", type, "\n")
  # }
  # Y = Y[Id]
  
  if (verbose) {
    cat("Proportion of positives after", type, ":", 
        round(length(which(Y == 1))/N * 100, digits = 2), "% of", N, "observations \n")
  }
  
  #transform the outout with the original labels
  Y <- factor(Y == 1, levels = c(FALSE, TRUE), labels = lev)
  
  return(list(X = X, Y = Y, id.rm = id.rm))
}
