#' @rdname rfPermute
#' 
#' @importFrom parallel makeForkCluster parLapply stopCluster
#' @export rfPermute.default
#' @export
#' 
rfPermute.default <- function(x, y, ..., nrep = 100, num.cores = 1) {  
  orig.call <- match.call()
  orig.call$nrep <- NULL
  orig.call$num.cores <- NULL
  orig.call[[1]] <- as.name("randomForest")
  rf.call <- orig.call
  rf.call$x <- x
  rf.call$y <- y
  imp.element <- pmatch(names(rf.call), "importance")
  if(!all(is.na(imp.element))) {
    imp.id <- which(!is.na(imp.element))
    names(rf.call)[imp.id] <- "importance"
  }    
  rf.call$importance <- TRUE
  rf.call[-1] <- lapply(as.list(rf.call[-1]), eval, envir = parent.frame())
  rf <- eval(rf.call)
  rf$call <- orig.call
  
  # permutes 'y' in rf.call 'nrep' times and runs randomForest  
  if(nrep > 0) {
    # define permutation function
    #   returns a 3-dimensional array of unscaled and scaled importance scores
    permFunc <- function(y, perm.rf.call) {
      rf.call$y <- y
      perm.rf <- eval(rf.call)
      imp <- perm.rf$importance
      impSD <- perm.rf$importanceSD
      rm(perm.rf)
      imp.arr <- .makeImpArray(imp, 2, c("unscaled", "scaled"))
      imp.arr[, , 2] <- .scaleImp(imp, impSD)
      return(imp.arr)
    }
    
    # create list of permuted y values
    ran.y <- lapply(1:nrep, function(i) sample(rf.call$y))
    
    # get importance scores for permutations
    #  a list of 3-dimensional arrays of importance scores
    null.dist<- NULL
    if(num.cores > 1) {
      cl <- makeForkCluster(num.cores)
      tryCatch({
        null.dist <- parLapply(cl, ran.y, permFunc, perm.rf.call = rf.call)
      }, finally = {
        stopCluster(cl)
        closeAllConnections()
      })
    } else {
      null.dist <- lapply(ran.y, permFunc)
    }
    
    # create and load null distribution arrays for each scaled and unscaled importances
    rf$null.dist <- list(unscaled = .makeImpArray(rf$importance, nrep, NULL))
    rf$null.dist$scaled <- rf$null.dist$unscaled
    for(i in 1:nrep) {
      rf$null.dist$unscaled[, , i] <- null.dist[[i]][, , "unscaled"]
      rf$null.dist$scaled[, , i] <- null.dist[[i]][, , "scaled"]
    }
    
    # calculate p-value of observed importance metrics
    rf$pval <- .calcImpPval(rf) 
    
    class(rf) <- c("rfPermute", "randomForest")
  }
  
  return(rf)  
}