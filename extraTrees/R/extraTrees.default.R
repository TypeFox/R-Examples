

## converts Java org.extratrees.data.Matrix to R matrix (doubles):
toRMatrix <- function( javam ) {
    return( matrix( 
        .jfield(javam, "[D", "v" ),
        .jfield(javam, "I", "nrows" ),
        .jfield(javam, "I", "ncols" )
    ))
}

## converts R matrix (or data.frame) into Java matrix (doubles)
toJavaMatrix <- function( m ) {
    if ( !is.matrix(m) ) { m = as.matrix(m) }
    if ( !is.double(m) ) { m = m + 0.0 }
    return(.jnew(
        "org.extratrees.data.Matrix", .jarray(m), nrow(m), ncol(m)
    ))
}

## converts R matrix (or data.frame) into Java CSparseMatrix
toJavaCSMatrix <- function( m ) {
    nzi = Matrix::which(m != 0, arr.ind=TRUE)
    v   = m[nzi]
    return(.jnew("org.extratrees.data.CSparseMatrix", 
                 .jarray(nzi[,1] - as.integer(1)),
                 .jarray(nzi[,2] - as.integer(1)),
                 .jarray(as.double(v)),
                 nrow(m),
                 ncol(m)
    ))
}

## converts R matrix, data.frame or sparseMatrix into Java Array2D
toJavaMatrix2D <- function( m ) {
  if ( is(m, "sparseMatrix") ) {
    jm = .jcast(toJavaCSMatrix(m), new.class="org/extratrees/data/Array2D")
  } else {
    jm = .jcast(toJavaMatrix(m), new.class="org/extratrees/data/Array2D")
  }
  return( jm )
}

## computes two integers, passed to Java for a long
get64BitSeed <- function() {
  as.integer(
    sample.int(n=2^31-1, size=2, replace=TRUE) * 
      c(-1, 1)[sample.int(n=2, size=2, replace=TRUE)]
  )
}

## creates a new extraTrees object based on selection
selectTrees <- function( object, selection ) {
    ## checking if right object:
    if (!inherits(object, "extraTrees")) {
        stop("object not of class extraTrees")
    }
    ## checking selection is logical:
    if (!is.logical(selection)) {
        stop("selection should be list of logical (T/F) values.")
    }
    ## checking selection has correct length:
    if (length(selection)!=object$ntree) {
        stop(sprintf("Length of selection (%d) should be equal to the number of trees (%d)", length(selection), object$ntree))
    }
    ## choosing the correct class:
    if (object$factor) {
        etClass = "Lorg/extratrees/FactorExtraTrees;"
    } else {
        etClass = "Lorg/extratrees/ExtraTrees;"
    }
    ## copying S3 object and creating new java object:
    etNew = object
    etNew$jobject = .jcall(object$jobject, etClass, "selectTrees", .jarray(selection) )
    etNew$ntree = .jcall(etNew$jobject, "I", "getNumTrees" )
    return(etNew)
}

prepareForSave <- function( object ) {
    if (!inherits(object, "extraTrees")) {
        stop("object not of class extraTrees")
    }
    .jcache(object$jobject)
}

## main extraTree training function
extraTrees.default <- function(x, y, 
             #xtest=NULL, ytest=NULL, 
             ntree=500,
             mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
             numRandomCuts = 1,
             evenCuts = FALSE,
             numThreads = 1,
             quantile = F,
             weights = NULL,
             subsetSizes = NULL,
             subsetGroups = NULL,
             tasks = NULL,
             probOfTaskCuts = mtry / ncol(x),
             numRandomTaskCuts = 1,
             na.action = "stop",
             ...) {
    args <- list(...)
    if (length(args) > 0) stop(sprintf("Illegal argument '%s' to extraTrees.", names(args)[1]))
  
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    ## making sure no NAs:
    if ( any(is.na(y)) ) stop("Output vector y contains NAs.")
    if ( !is.null(tasks) && any(is.na(tasks)) ) stop("Task vector contains NAs.")
    
    if ( numThreads < 1 ) stop("numThreads has to be 1 or bigger.")

    ## uncomment when xtest/ytest are used:
    #testdat <- !is.null(xtest)
    #if (testdat) {
    #    if (ncol(x) != ncol(xtest))
    #        stop("x and xtest must have same number of columns")
    #    ntest <- nrow(xtest)
    #    xts.row.names <- rownames(xtest)
    #}
    
    et <- list()
    et$ntree    = ntree
    et$nodesize = nodesize
    et$ndim     = p
    et$mtry     = mtry
    et$factor   = is.factor(y)
    et$numRandomCuts = numRandomCuts
    et$evenCuts = evenCuts
    et$numThreads = numThreads
    et$quantile   = quantile
    et$useWeights = ! is.null(weights)
    et$usesubsetting = ! is.null(subsetSizes) && sum(subsetSizes) != nrow(x)
    et$multitask  = ! is.null(tasks)
    et$probOfTaskCuts = probOfTaskCuts
    et$numRandomTaskCuts = numRandomTaskCuts
    et$call <- match.call()
    et$call[[1]] <- as.name("extraTrees")
    
    class(et) = "extraTrees"

    if (nrow(x) != length(y)) {
        stop(sprintf("Length of y (%d) is not equal to the number of samples in x (%d).", length(y), nrow(x) ) )
    }
    
    if ( ! et$factor && ! is.numeric(y) ) {
      stop("y values have to be either factor (classification) or numeric (regression).")
    }

    et$xHasNA = FALSE
    if (na.action == "stop") {
      if ( any(is.na(x)) ) stop("Input matrix x contains NAs. Change na.action to 'zero' or 'fuse' to allow NA to be used or remove samples with NA.")
    } else if (na.action == "zero") {
      x[ is.na(x) ] = 0
    } else if (na.action == "fuse") {
      et$xHasNA = any(is.na(x))
    } else {
      stop("na.action should be either 'stop', 'zero' or 'fuse'. See manual for details.")
    }
    
    if ( ! is.null(weights) ) {
      if (nrow(x) != length(weights)) {
        stop(sprintf("Length of weights (%d) is not equal to the number of samples in x (%d).", length(weights), nrow(x) ) )
      }
      if (any(weights <= 0)) {
        stop("Weights have to be positive.")
      }
    }
    
    if ( et$usesubsetting ) {
      if (sum(subsetSizes) > nrow(x)) {
        stop(sprintf("Total of subsetSizes (%d) should not be bigger than the number of samples in x (%d).", sum(subsetSizes), nrow(x) ))
      }
      ## if only one subset size then subsetGroups are not used
      if (length(subsetSizes) >= 2) {
        if (nrow(x) != length(subsetGroups)) {
          stop(sprintf("Length of subsetGroups (%d) is not equal to the number of samples in x (%d).", length(weights), nrow(x) ) )
        }
        if ( ! is.factor(subsetGroups) ) {
          subsetGroups = as.factor( subsetGroups )
        }
        
        numUnique = length(levels(subsetGroups))
        if (numUnique != length(subsetSizes)) {
          stop(sprintf("Number of unique subsetGroups (%d) has to be the same as length of number of subsetSizes (%d).", numUnique, length(subsetSizes)  ))
        }
      }
    }
    
    ## making sure if tasks is present there are only two factors
    if ( ! is.null(tasks) ) {
        if (nrow(x)!=length(tasks)) {
            stop(sprintf("Length of tasks (%d) is not equal to the number of inputs in x (%d).", length(tasks), nrow(x) ) )
        }
        if ( et$factor && length(unique(y)) != 2 ) {
            stop("Multi-task learning only works with 2 factors (binary classification). 3 or more classes is not supported.")
        }
        if (min(tasks) < 1) {
            stop("Tasks should be positive integers.")
        }
        if ( et$quantile ) {
            stop("Quantile regression is not (yet) supported with multi-task learning.")
        }
        if ( et$useWeights ) {
            stop("Weights are not (yet) supported with multi-task learning.")
        }
    }

    if (et$factor && length(unique(y)) < 2) {
        stop("Need at least two classes to do classification.")
    }
    if (et$factor) {
        if (quantile) {
            stop("Option quantile cannot be used for classification.")
        }
        ## classification:
        et$levels = levels(y)
        ## creating FactorExtraTree object with the data
        et$jobject = .jnew(
            "org.extratrees.FactorExtraTrees",
            toJavaMatrix2D(x),
            .jarray( as.integer( as.integer(y)-1 ) )
        )
        .jcall( et$jobject, "V", "setnFactors", as.integer(length(et$levels)) )
    } else if (et$quantile) {
        ## quantile regression:
        et$jobject = .jnew(
            "org.extratrees.QuantileExtraTrees",
            toJavaMatrix2D(x),
            .jarray( as.double(y) )
        )
    } else {
        ## regression:
        ## creating ExtraTree object with the data
        et$jobject = .jnew(
            "org.extratrees.ExtraTrees",
            toJavaMatrix2D(x),
            .jarray( as.double(y) )
        )
    }
    
    #if (et$xHasNA) {
    #  print("Note: Input matrix has NA. See ?extraTrees how ET builds trees with NAs.")
    #}
    
    ## setting variables:
    .jcall( et$jobject, "V", "setNumRandomCuts", as.integer(et$numRandomCuts) )
    .jcall( et$jobject, "V", "setEvenCuts", et$evenCuts )
    .jcall( et$jobject, "V", "setNumThreads", as.integer(et$numThreads) )
    .jcall( et$jobject, "V", "setHasNaN", et$xHasNA)
    seeds = get64BitSeed()
    .jcall( et$jobject, "V", "setSeed", seeds[1], seeds[2])
    
    ## if present set weights:
    if (et$useWeights) {
      .jcall( et$jobject, "V", "setWeights", .jarray(as.double(weights)) )
    }
    
    ## if given set subsetting:
    if (et$usesubsetting) {
      if (length(subsetSizes) == 1) {
        .jcall( et$jobject, "V", "setSubsetting", as.integer(subsetSizes[1]) )
      } else {
        .jcall( et$jobject, "V", "setSubsetting", 
                .jarray(as.integer( subsetSizes )), 
                .jarray(as.integer( as.integer(subsetGroups)-1 )) 
              )
      }
    }
    
    ## multitask variables:
    if (et$multitask) {
        .jcall( et$jobject, "V", "setTasks", .jarray(as.integer(tasks-1)) )
        .jcall( et$jobject, "V", "setProbOfTaskCuts", et$probOfTaskCuts )
        .jcall( et$jobject, "V", "setNumRandomTaskCuts", as.integer(et$numRandomTaskCuts) )
    }
    
    ## learning the trees (stored at the et$jobject)
    .jcall( et$jobject, "V", "learnTrees", as.integer(et$nodesize), as.integer(et$mtry), as.integer(et$ntree) )
    
    return( et )
}

