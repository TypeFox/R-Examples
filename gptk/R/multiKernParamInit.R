multiKernParamInit <-
function (kern) {

  kern$nParams <- 0
  kern$transforms <- list()

  if ( !("comp" %in% names(kern)) )
    kern$comp <- list()

  kern$numBlocks <- length(kern$comp)
  kern$isStationary <- TRUE

  kern$block <- list()
  for ( i in seq(along=kern$comp) ) {
    if ( !kern$comp[[i]]$isStationary )
      kern$isStationary <- FALSE

    kern$comp[[i]] <- kernParamInit(kern$comp[[i]])
    kern$nParams <- kern$nParams + kern$comp[[i]]$nParams
    kern$comp[[i]]$index <- array()

    kern$block[[i]] <- list(cross=array(), transpose=array())

    for ( j in seq(length.out=i-1) ) {
      if ( .kernTestCombinationFunction(kern$comp[[i]], kern$comp[[j]]) ) {
        kern$block[[i]]$cross[j] <- paste(kern$comp[[i]]$type, "X", kern$comp[[j]]$type, sep="")
        kern$block[[i]]$transpose[j] <- FALSE
      } else {
        if ( .kernTestCombinationFunction(kern$comp[[j]], kern$comp[[i]]) ) {
          kern$block[[i]]$cross[j] <- paste(kern$comp[[j]]$type, "X", kern$comp[[i]]$type, sep="")
          kern$block[[i]]$transpose[j] <- TRUE
        } else {
          warning(paste("No cross covariance found between", kern$comp[[i]]$type, "and", kern$comp[[j]]$type, "assuming independence."))
          kern$block[[i]]$cross[j] <- ""
          kern$block[[i]]$transpose[j] <- 0
        }
      }
    }
  }

  kern$paramGroups <- diag(1, nrow=kern$nParams, ncol=kern$nParams)

  kern$fixBlocks <- rep(FALSE, kern$numBlocks)
  if ("options" %in% names(kern) && "fixedBlocks" %in% names(kern$options)
      && kern$options$fixedBlocks) {
    kern$fixBlocks[kern$options$fixedBlocks] <- TRUE
    kern$cache <- new.env(parent=emptyenv())
    assign("cache", list(list(list())), envir=kern$cache)
  }

  return (kern)
}
