#'@title Monte Carlo simulations with matrixpls
#'
#'@description
#'Performs Monte Carlo simulations of \code{\link{matrixpls}} with the \code{\link[simsem]{sim}} function of the \code{simsem} package.
#'The standard errors and confidence intervals are estimated with the \code{\link[boot]{boot}} and \code{\link[boot]{boot.ci}} functions
#'of the \code{boot} package.
#
#'@details
#' This funtion calls the \code{\link[simsem]{sim}} function from the \code{simsem} package to perform Monte
#' Carlo simulations with matrixpls. The function parses the model parameters and replaces it with
#' a function call that estimates the model and bootstrapped standard errors and confidence
#' intervals with \link{matrixpls.boot}.
#' 
#' If the \code{generate} or \code{rawdata} arguments are not specified in the \code{\link[simsem]{sim}} arguments
#'  then the \code{model} argument will be used for data generation and must be specified in lavaan format.
#'
#'@param nRep Number of replications. If any of the \code{n}, \code{pmMCAR}, or \code{pmMAR} arguments are specified as lists, the number of replications will default to the length of the list(s), and \code{nRep} need not be specified.
#'
#'@inheritParams matrixpls
#'
#'@param ... All other arguments are passed through to \code{\link[simsem]{sim}},
#' \code{\link{matrixpls.boot}}, or \code{\link{matrixpls}}.
#'
#'@param cilevel Confidence level. This argument will be forwarded to the \code{\link[boot]{boot.ci}} when calculating the confidence intervals.
#'
#'@param citype Type of confidence interval. This argument will be forwarded to the \code{\link[boot]{boot.ci}} when calculating the confidence intervals.
#'
#'@param boot.R Number of bootstrap replications to use to estimate standard errors or \code{FALSE} to disable bootstrapping.
#'
#'@param fitIndices A function that returns a list of fit indices for the model. Setting this argument to \code{NULL} disables fit indices.
#'
#'@param prefun A function to be applied to the dataset before each replication. The output of this
#'function is passed as arguments to \code{\link{matrixpls}}
#'
#'@param outfun A function to be applied to the  \code{matrixpls} output at each replication. 
#'Output from this function in each replication will be saved in the simulation 
#'output (SimResult), and can be obtained using the getExtraOutput function.

#'@param outfundata A function to be applied to the \code{matrixpls} output and the 
#'generated data after each replication. Users can get the characteristics of the 
#'generated data and also compare the characteristics with the generated output. 
#'The output from this function in each replication will be saved in the 
#'simulation output (SimResult), and can be obtained using the getExtraOutput function.
#'
#'@inheritParams simsem::sim
#'
#'@return An object of class \code{\link[simsem]{SimResult-class}}.
#'
#'@example example/fragment-requireSimsem.R
#'@example example/matrixpls.sim-example2.R
#'@example example/fragment-endBlock.R
#'
#'@seealso
#'
#'\code{\link{matrixpls}}, \code{\link{matrixpls.boot}}, \code{\link[simsem]{sim}}, \code{\link[simsem]{SimResult-class}}
#'
#'@include matrixpls.R
#'
#'@export


matrixpls.sim <- function(nRep = NULL, model = NULL, n = NULL, ..., cilevel = 0.95,
                          citype=c("norm","basic", "stud", "perc", "bca"), 
                          boot.R = 500, fitIndices = fitSummary,
                          outfundata = NULL, outfun = NULL,
                          prefun = NULL){
  
  if(! requireNamespace("simsem")) stop("matrixpls.sim requires the simsem package")
  
  
  # Basic verification of the arguments
  if(boot.R != FALSE)	assertive::assert_all_are_positive(boot.R)	
  assertive::assert_all_are_positive(cilevel)
  assertive::assert_all_are_true(cilevel<1)
  citype <- match.arg(citype)
  
  # Decide which arguments to route to simsem and which to matrispls.boot
  
  allArgs <- list(...)	
  simsem.argNames <- names(formals(simsem::sim))
  
  whichArgsToSimsem <- names(allArgs) %in% simsem.argNames
  simsemArgs <- allArgs[whichArgsToSimsem]
  matrixplsArgs <- allArgs[! whichArgsToSimsem]
  
  nativeModel <- parseModelToNativeFormat(model)
  
  # Because we are using a custom estimator function, the generate model must be set.
  
  if(!"generate" %in% names(simsemArgs) &&
       !"rawData" %in% names(simsemArgs)){
    
    # Derive a SimSem model from lavaan parameter table by fitting a model to a diagonal matrix
    
#    cvMat <- diag(nrow(nativeModel$reflective))
#    colnames(cvMat) <- rownames(cvMat) <- rownames(nativeModel$reflective)
#    fit <- lavaan::lavaan(model, sample.cov = cvMat, sample.nobs = 100)
#    simsemArgs$generate <- simsem::model.lavaan(fit)
    simsemArgs$generate <- model
  }
  
  if(!"W.mod" %in% names(matrixplsArgs)) matrixplsArgs$W.mod<- defaultWeightModelWithModel(model)
  
  matrixplsArgs <- c(list(R= boot.R, model = nativeModel),matrixplsArgs)
  
  #
  # The LV scores returned by SimSem for endogenous LVs are not the LVs, but their error terms.
  # To get scores for the endogenous LVs, we need to calculate these based on the other LVs. To
  # so this, the population model needs to be parsed.
  #
  
  partable <- NULL
  
  if(is.character(model)) {
    # Remove all multigroup specifications because we do not support multigroup analyses
    model <- gsub("c\\(.+?\\)","NA",model)
    partable <- lavaan::lavaanify(model)
  } else if (is.partable(model)) {
    partable <- model
  } else if (methods::is(model, "lavaan")) {
    partable <- model@ParTable
  }
  
  if(!is.null(partable)){
    
    factorLoadings <- partable[partable$op == "=~",]
    regressions <- partable[partable$op == "~",]
    formativeLoadings <- partable[partable$op == "<~",]
    
    # Parse the variables
    latentVariableNames <- unique(c(factorLoadings$lhs, formativeLoadings$lhs))
    
    # Set up empty model tables
    
    inner <- matrix(0,length(latentVariableNames),length(latentVariableNames))
    colnames(inner)<-rownames(inner)<-latentVariableNames
    
    # Set the relationships in the tables
    
    latentRegressions <- regressions[regressions$rhs %in% latentVariableNames & 
                                       regressions$lhs %in% latentVariableNames,]
    
    rows <- match(regressions$lhs, latentVariableNames)
    cols <- match(regressions$rhs, latentVariableNames)
    indices <- rows + (cols-1)*nrow(inner) 
    
    inner[indices] <- regressions$ustart
    
  }
  
  # A function that takes a data set and returns a list. The list must
  # contain at least three objects: a vector of parameter estimates (coef),
  # a vector of standard error (se), and the convergence status as TRUE or
  # FALSE (converged). There are five optional objects in the list: a vector
  # of fit indices (fit), a vector of standardized estimates (std), any
  # extra output (extra), fraction missing type I (FMI1), and fraction
  # missing type II (FMI2). Note that the coef, se, std, FMI1, and FMI2 must
  # be a vector with names. The name of those vectors across different
  # objects must be the same.
  
  
  modelFun  <- function(data){

    # Indices for parameters excluding weights
    
    parameterIndices <- 1:(sum(nativeModel$inner) + sum(nativeModel$reflective) + sum(nativeModel$formative))

    
    # Convert the data to matrix for efficiency
    if(boot.R == FALSE){
      S <- stats::cov(data)
      
      if(! is.null(prefun)){
        extraArgs <- prefun(data)
      }
      else extraArgs <- list()
      
      matrixpls.res <- do.call(matrixpls, c(list(S), matrixplsArgs,extraArgs))
    }
    else{
      boot.out <- do.call(matrixpls.boot, c(list(as.matrix(data)), matrixplsArgs, list(prefun = prefun)))
      matrixpls.res  <- boot.out$t0
    }
    
    # Check for inadmissible solutions. 
    #
    # 0: Converged normally 
    # 1: Non-convergent result
    # 2: Non-converged imputation (not used)
    # 3: At least one SE is negative or NA (not used)
    # 4: At least one variance estimate is negative
    # 5: At least one correlation estimate is greater than 1 or less than -1
    
    # Non-iterative weight functiosn do not return convergence status so both NULL
    # and TRUE are considered as converged
    
    if(is.null(attr(matrixpls.res,"converged")) ||
       attr(matrixpls.res,"converged")){
      
      converged <- 0
      
      C <- attr(matrixpls.res,"C")
      
      if(max(abs(C[lower.tri(C)]))>1){
        converged <- 5 
      }
      # If the model is estimated with 2SLS, then checking C is not enough to check for admissible
      # solution. We need to calculate the explained variances of the endogenous composites
      
      else{
        inner <- attr(matrixpls.res,"inner")
        if(any(diag(inner%*%C%*%t(inner)) > 1)){
          converged <- 4
        }
      }
    }
    else converged <- 1
    
    ret <- list(coef = matrixpls.res[parameterIndices],
                converged = converged)
    
    # If the data were generated sequentially using LV scores, calculate the true reliabilities
    
    latentVar <- attr(data,"latentVar")
    
    if(! is.null(latentVar)){
      lvScores <-  as.matrix(data) %*% t(attr(matrixpls.res, "W"))
      
      # The latent vars and composites should be in the same order. 
      trueScores <- as.matrix(latentVar[,1:ncol(lvScores)])
      
      r <- diag(stats::cor(lvScores,trueScores))
      
      # Keep the sign of the correlation when calculating reliabilities
      R <- sign(r) * r^2
      names(R) <- colnames(lvScores)
      attr(matrixpls.res, "R") <- R
      
    }
    
    # Store CIs and SEs if bootstrapping was done
    
    if(boot.R != FALSE){
      
      cis <- sapply(parameterIndices, FUN = function(index) {
        boot.ci.out <- boot::boot.ci(boot.out, conf = cilevel, type = citype, index=index)
        
        # The cis start from the fourth slot and we only have one type of ci. 
        # The list names do not match the type parameter exactly (e.g. "norm" vs. "normal")
        
        cis <- boot.ci.out[[4]]
        cis[,ncol(cis)-1:0]
      })
      
      ses <- apply(boot.out$t[,parameterIndices],2,stats::sd)
      names(ses) <- names(ret[["coef"]])
      colnames(cis) <- names(ret[["coef"]])
      
      ret <- c(ret, list(se = ses,
                         cilower = cis[1,],
                         ciupper = cis[2,],
                         extra = boot.out))
      
    }
    else{
      # simsem requires some se estimates, so return a matrix of NAs
      ses <- rep(NA, length(parameterIndices))
      names(ses) <- names(ret[["coef"]])
      ret <- c(ret, list(se = ses,
                         extra = matrixpls.res))
      
    }
    if(! is.null(fitIndices)){
      assertive::assert_is_function(fitIndices)
      fitlist <- unlist(fitIndices(matrixpls.res))
      ret$fit <- fitlist
    }
    else ret$fit <- c()
    
    # Apply the outfun
    if(!is.null(outfundata)){
      
      if(! is.null(outfun)) warning("Both outfundata and outfun were specified. Only outfundata is applied")
      ret$extra <- outfundata(ret$extra, data)
    }
    else if(!is.null(outfun)){
      ret$extra <- outfun(ret$extra)      
    }
    
    return(ret)
  }
  
  simsemArgs <- c(list(nRep = nRep, model = modelFun, n = n,
                       outfun = function(x){x$extra}),
                  simsemArgs)

  ret <- do.call(simsem::sim, simsemArgs)
  
  ret
}
