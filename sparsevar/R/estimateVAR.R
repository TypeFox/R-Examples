#' @title Multivariate VAR estimation
#' 
#' @description A function to estimate a (possibly high-dimensional) multivariate VAR time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#' @param data the data from the time series: variables in columns and observations in 
#' rows
#' @param p order of the VAR model
#' @param penalty the penalty function to use. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}
#' @param options list of options for the function. Global options are:
#' \code{threshold}: \code{TRUE} or \code{FALSE}. If \code{TRUE} all the entries smaller 
#' than the oracle threshold are set to zero. \code{scale} scale the data?
#' 
#' @return A the list (of length \code{p}) of the estimated matrices of the process
#' @return fit the results of the penalized LS estimation
#' @return mse the mean square error of the cross validation
#' @return time elapsed time for the estimation
#'
#' @usage estimateVAR(data, p = 1, penalty = "ENET", options = NULL)
#' 
#' @export
estimateVAR <- function(data, p = 1, penalty = "ENET", options = NULL) {

  # get the number of rows and columns
  nr <- nrow(data)
  nc <- ncol(data)
  
  # make sure the data is in matrix format
  data <- as.matrix(data)

  # scale the matrix columns
  scale <- ifelse(is.null(options$scale), FALSE, options$scale)  
  if (scale == TRUE) {
    data <- apply(FUN = scale, X = data, MARGIN = 2)
  }
  
  # create Xs and Ys (temp variables)
  tmpX <- data[1:(nr-1), ]
  tmpY <- data[2:(nr), ]
  
  # create the data matrix
  tmpX <- duplicateMatrix(tmpX, p)
  tmpY <- tmpY[p:nrow(tmpY), ]

  y <- as.vector(tmpY)
  
  # Hadamard product for data
  I <- Matrix::Diagonal(nc)
  X <- kronecker(I, tmpX)
  
  if (penalty == "ENET") {
    
    # By default repeatedCV = FALSE
    options$repeatedCV <- ifelse(is.null(options$repeatedCV), FALSE, TRUE)
    
    # fit the ENET model
    t <- Sys.time()
    fit <- varENET(X, y, options)
    elapsed <- Sys.time() - t
    
    if (options$repeatedCV == FALSE){
      # extract what is needed
      lambda <- ifelse(is.null(options$lambda), "lambda.min", options$lambda)
      # extract the coefficients and reshape the matrix
      Avector <- stats::coef(fit, s = lambda)
      A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    } else {
      Avector <- fit$Avector
      A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    }
    
    mse <- min(fit$cvm)
    
  } else if (penalty == "SCAD") {
    
    # convert from sparse matrix to std matrix (SCAD does not work with sparse matrices)
    X <- as.matrix(X)
    # fit the SCAD model
    t <- Sys.time()
    fit <- varSCAD(X, y, options)
    elapsed <- Sys.time() - t
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    mse <- min(fit$cve)
    
  } else if (penalty == "MCP") {
    
    # convert from sparse matrix to std matrix (MCP does not work with sparse matrices)
    X <- as.matrix(X)
    # fit the MCP model
    t <- Sys.time()
    fit <- varMCP(X, y, options)
    elapsed <- Sys.time() - t
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    mse <- min(fit$cve)
    
  } else {
    
    # Unknown penalty error
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
    
  }
  
  # If threshold = TRUE then set to zero all the entries that are small
  if (!is.null(options$threshold)) {
    if (options$threshold == TRUE) {
      tr <- 1 / sqrt(p*nc*log(nr))
      L <- abs(A) >= tr
      A <- A * L
    }
  }
  
  # Get back the list of VAR matrices (of length p)
  A <- splitMatrix(A, p)
  
  # Output
  output = list()
  output$A <- A
  output$fit <- fit
  output$mse <- mse
  output$time <- elapsed
  return(output)
  
}

varENET <- function(X,y, options = NULL) {
  
  a  <- ifelse(is.null(options$alpha), 1, options$alpha)
  nl <- ifelse(is.null(options$nlambda), 100, options$nlambda)
  tm <- ifelse(is.null(options$type.measure), "mse", options$type.measure)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)
  repeatedCV <- options$repeatedCV #ifelse(is.null(options$repeatedCV), FALSE, TRUE)
  nRepeats <- ifelse(is.null(options$nRepeats), 3, options$nRepeats)
  
  if (repeatedCV == TRUE) {
    
    trCtrl <- caret::trainControl(method = "repeatedcv", number = nf, repeats = nRepeats)
    lam <- glmnet::glmnet(X, y, alpha = a)$lambda
    gr <- expand.grid(.alpha = a, .lambda = lam)
    fit <- caret::train(x = X, y = y, method = "glmnet", trControl = trCtrl, tuneGrid = gr)
    b <- stats::coef(fit$finalModel, fit$bestTune$lambda)
    cvm <- 3
    fit <- list()
    fit$Avector <- b
    fit$cvm <- cvm
    return(fit)
    
  } 
  
  # Assign ids to the CV-folds (useful for replication of results)  
  if (is.null(options$foldsIDs)) {
    foldsIDs <- numeric(0)
  } else {
    nr <- nrow(X)
    foldsIDs <- sort(rep(seq(nf), length = nr))
  }
  
  ##############################################################################
  # TODO: change parallel backend (parallel -> doMC)
  ##############################################################################
  if(parall == TRUE) {
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      #cl <- registerDoMC(ncores)
      cl <- parallel::makeCluster(ncores)
      if (length(foldsIDs) == 0) {
        cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = TRUE)
      } else {
        cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, foldid = foldsIDs, parallel = TRUE)
      }
      parallel::stopCluster(cl)
    }
  } else {
    if (length(foldsIDs) == 0) {
      cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = FALSE)
    } else {
      cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, foldid = foldsIDs, parallel = FALSE)
    }
    
  }
  
  return(cvfit)
  
}

varSCAD <- function(X, y, options = NULL) {
  
  e <- ifelse(is.null(options$eps), 0.01, options$eps)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)
  
  if(parall == TRUE) {
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e, cluster = cl)
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e)
  }

  return(cvfit)
  
}

varMCP <- function(X, y, options = NULL) {
  
  e <- ifelse(is.null(options$eps), 0.01, options$eps)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)
  
  if(parall == TRUE) {
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e, cluster = cl)
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e)
  }
  
  return(cvfit)
  
}


splitMatrix <- function(M, p) {
  
  nr <- nrow(M)
  A <- list()
  
  for (i in 1:p) {
    
    ix <- ((i-1) * nr) + (1:nr)
    A[[i]] <- M[1:nr, ix]  
    
  }

  return(A)
}

duplicateMatrix <- function(data, p) {
  
  nr <- nrow(data)
  nc <- ncol(data)
  
  outputData <- data
  
  if (p > 1) {
    for (i in 1:(p-1)) {
      
      tmpData <- matrix(0, nrow = nr, ncol = nc)
      tmpData[(i+1):nr, ] <- data[1:(nr-i), ]
      outputData <- cbind(outputData, tmpData)
      
    }
  }
  
  outputData <- outputData[p:nr, ]
  return(outputData)
  
}
