# =========== Main functions ===========

#'@title A semPLS compatibility wrapper for matrixpls
#'
#'@description
#'\code{matrixpls.sempls} mimics \code{\link[semPLS]{sempls}} function of the \code{semPLS} package.
#'The arguments and their default values and the output of the function are identical with \code{\link[semPLS]{sempls}} function,
#'but internally the function uses matrixpls estimation.
#'
#'@param model An object inheriting from class \code{plsm} as returned from \code{\link[semPLS]{plsm}}
#' or \code{\link[semPLS]{read.splsm}}.  
#' 
#'@param ... Other arguments are ignored
#'
#'@inheritParams semPLS::sempls
#'
#'@return An object of class \code{\link[semPLS]{sempls}}. 
#'
#'@references Monecke, A., & Leisch, F. (2012). semPLS: Structural Equation Modeling Using Partial Least Squares. \emph{Journal of Statistical Software}, 48(3), 1–32.

#'
#'@seealso
#'\code{\link[semPLS]{sempls}}
#'
#'@export
#'@example example/fragment-requireSemPLS.R
#'@example example/matrixpls.sempls-example.R
#'@example example/fragment-endBlock.R


matrixpls.sempls <-
  function(model, data, maxit=20, tol=1e-7, scaled=TRUE, sum1=FALSE, wscheme="centroid", pairwise=FALSE,
           method=c("pearson", "kendall", "spearman"),
           convCrit=c("relative", "square"), verbose=TRUE, ...){
    
    # This function is largely copied from the semPLS package (licensed under GPL 2)
    # The code specific to matrixpls are marked in the following way
        
    method <- match.arg(method)
    convCrit <- match.arg(convCrit)
    result <- list(coefficients=NULL, path_coefficients=NULL,
                   outer_loadings=NULL ,cross_loadings=NULL,
                   total_effects=NULL,inner_weights=NULL, outer_weights=NULL,
                   blocks=NULL, factor_scores=NULL, data=NULL, scaled=scaled,
                   model=model, weighting_scheme=NULL, weights_evolution=NULL,
                   sum1=sum1, pairwise=pairwise, method=method, iterations=NULL,
                   convCrit=convCrit, verbose=verbose, tolerance=tol, maxit=maxit, N=NULL,
                   incomplete=NULL, Hanafi=NULL)
    class(result) <- "sempls"
    
    # checking the data
    data <- data[, model$manifest]
    N <- nrow(data)
    missings <- which(stats::complete.cases(data)==FALSE)
    if(length(missings)==0 & verbose){
      cat("All", N ,"observations are valid.\n")
      if(pairwise){
        pairwise <- FALSE
        cat("Argument 'pairwise' is reset to FALSE.\n")
      }
    }
    else if(length(missings)!=0 & !pairwise & verbose){
      # Just keeping the observations, that are complete.
      data <- stats::na.omit(data[, model$manifest])
      cat("Data rows:", paste(missings, collapse=", "),
          "\nare not taken into acount, due to missings in the manifest variables.\n",
          "Total number of stats::complete.cases:", N-length(missings), "\n")
    }
    else if(verbose){
      cat("Data rows", paste(missings, collapse=", "),
          " contain missing values.\n",
          "Total number of stats::complete.cases:", N-length(missings), "\n")
    }
    ## check the variances of the data
    if(!all(apply(data, 2, stats::sd, na.rm=TRUE) != 0)){
      stop("The MVs: ",
           paste(colnames(data)[which(apply(data, 2, stats::sd)==0)], collapse=", "),
           "\n  have standard deviation equal to 0.\n",
           "  Recheck model!\n")
    }
    
    ## scale data?
    # Note: scale() changes class(data) to 'matrix'
    if(scaled) data <- scale(data)
    
    #############################################
    # Select the function according to the weighting scheme
    if(wscheme %in% c("A", "centroid")) {
      # Start of matrixpls code
      innerEstimator <- inner.centroid
      # End of matrixpls code
      result$weighting_scheme <- "centroid"
    }
    else if(wscheme %in% c("B", "factorial")) {
      # Start of matrixpls code
      innerEstimator <- inner.factor
      # End of matrixpls code
      result$weighting_scheme <- "factorial"
    }
    else if(wscheme %in% c("C", "pw", "pathWeighting")) {
      # Start of matrixpls code
      innerEstimator <- inner.path
      # End of matrixpls code
      result$weighting_scheme <- "path weighting"
    }
    else {stop("The argument E can only take the values 'A', 'B' or 'C'.\n See ?sempls")}
    
    # Start of matrixpls code
    
    modes <- unlist(lapply(model$blocks, function(x) attr(x,"mode")))
    modeA <- modes == "A"
    
    if(max(modeA) == 0) outerEstimators <- outer.modeB	
    else if(min(modeA) == 1) outerEstimators <- outer.modeA
    else{
      outerEstimators <- list(rep(NA,length(modeA)))
      outerEstimators[!modeA] <- list(outer.modeB)
      outerEstimators[modeA] <- list(outer.modeA)
    }
    
    S <- stats::cov(data)
    
    if(convCrit=="relative"){
      convCheck <- convCheck.relative
    }
    else if(convCrit=="square"){
      convCheck <- convCheck.square
    }	
    
    
    matrixpls.model <- list(inner = t(model$D), 
                            reflective = model$M, 
                            formative = matrix(0, ncol(model$M), nrow(model$M),
                                               dimnames = list(colnames(model$M), rownames(model$M))))
    
    matrixpls.res <- matrixpls(S,
                               matrixpls.model,
                               innerEstimator = innerEstimator,
                               outerEstimators = outerEstimators,
                               convCheck = convCheck, tol = tol,
                               standardize = FALSE)
    
    converged <- attr(matrixpls.res, "converged")
    i <- attr(matrixpls.res, "iterations") + 1
    
    # Create a new matrix to clear attribute
    W <- t(attr(matrixpls.res, "W"))
    Wnew <- matrix(W,nrow(W), ncol(W))
    colnames(Wnew) <- colnames(W)
    rownames(Wnew) <- rownames(W)
    
    innerWeights <- t(attr(matrixpls.res, "E"))
    whist <- attr(matrixpls.res, "history")
    
    weights_evolution <- apply(whist, 1, function(x){
      
      Wnew[Wnew != 0] <- x
      ret <- stats::reshape(as.data.frame(Wnew),
                     v.names="weights",
                     ids=rownames(Wnew),
                     idvar="MVs",
                     times=colnames(Wnew),
                     timevar="LVs",
                     varying=list(colnames(Wnew)),
                     direction="long")
      ret
    })
    
    weights_evolution <- do.call(rbind,weights_evolution)
    iterationCol <- data.frame(iteration = rep(0:(nrow(whist)-1), each = length(Wnew)))
    weights_evolution <- cbind(weights_evolution, iterationCol)
    
    # End of matrixpls code
    
    
    ## print
    if(converged & verbose){
      cat(paste("Converged after ", (i-1), " iterations.\n",
                "Tolerance: ", tol ,"\n", sep=""))
      if (wscheme %in% c("A", "centroid")) cat("Scheme: centroid\n")
      if (wscheme %in% c("B", "factorial")) cat("Scheme: factorial\n")
      if (wscheme %in% c("C", "pw", "pathWeighting")) cat("Scheme: path weighting\n")
    }
    else if(!converged){
      stop("Result did not converge after ", result$maxit, " iterations.\n",
           "\nIncrease 'maxit' and rerun.", sep="")
    }
    
    weights_evolution <- weights_evolution[weights_evolution!=0,]
    weights_evolution$LVs <- factor(weights_evolution$LVs,  levels=model$latent)
    # create result list
    ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
    
    # Start of matrixpls code
    result$path_coefficients <- t(attr(matrixpls.res,"inner"))
    result$cross_loadings <- t(attr(matrixpls.res,"IC"))
    # End of matrixpls code
    
    result$outer_loadings <- result$cross_loadings
    result$outer_loadings[Wnew==0] <- 0
    
    # Start of matrixpls code
    result$total_effects <- result$path_coefficients
    endo <- which(apply(matrixpls.model$inner != 0,1,any))
    result$total_effects[,endo] <- t(effects(matrixpls.res)$Total)
    # End of matrixpls code
    
    result$inner_weights <- innerWeights
    result$outer_weights <- Wnew
    result$weights_evolution <- weights_evolution
    
    # Start of matrixpls code
    result$Hanafi <- 		apply(whist,1,function(x){
      
      W[W != 0] <- x
      C <- t(W) %*% S %*% W
      
      matrix(c(sum(abs(C) * model$D),
               sum(C^2 * model$D)), nrow = 1)
      
    })
    
    result$Hanafi <- t(result$Hanafi)
    result$Hanafi <- cbind(result$Hanafi, 0:(i-1))
    
    colnames(result$Hanafi) <- c("f", "g", "iteration")
    
    #
    # semPLS calculates the weights using non-standardized weights and standardizes afterwards.
    # because of this the attributes scaled:center and scaled:scale added by scale do not match
    # the semPLS output. These are of little interest to the end user.
    #
    
    result$factor_scores <- scale(data %*% W)
    
    # End of matrixpls code
    
    result$data <- data
    result$N <- N
    result$incomplete <- missings
    result$iterations <- (i-1)
    
    ref <- matrixpls.model$reflective
    inn <- t(matrixpls.model$inner)
    
    pathNames <- paste(colnames(ref)[col(ref)[ref==1]], "->", rownames(ref)[row(ref)[ref==1]])
    estimates <- t(attr(matrixpls.res,"IC"))[ref==1]
    coefNames <- paste("lam",col(ref)[ref==1],unlist(apply(ref,2,function(x) 1:sum(x))),sep="_")
    
    pathNames <- c(pathNames, paste(rownames(inn)[row(inn)[inn==1]], "->", colnames(inn)[col(inn)[inn==1]]))
    estimates <- c(estimates, t(attr(matrixpls.res,"inner"))[inn==1])
    coefNames <- c(coefNames, paste("inner",row(inn)[inn==1],col(inn)[inn==1],sep="_"))
    
    # Start of matrixpls code
    
    result$coefficients <- data.frame(Path = pathNames,
                                      Estimate = estimates, 
                                      row.names = coefNames)
    # End of matrixpls code
    
    return(result)
  }


