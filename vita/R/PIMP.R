PIMP = function(X, ...)  UseMethod("PIMP")

PIMP.default = function(X, y, rForest, S = 100, parallel = FALSE, ncores=0, seed = 123, ...){

  ##################################################################
  # randomForest?
  if (!inherits(rForest, "randomForest")) stop("rForest is not of class randomForest")
  mtry=rForest$mtry
  ntree=rForest$ntree
  n = nrow(X)
  p = ncol(X)

  ##################################################################
  # Windows and parallel ?
  if(Sys.info()[['sysname']] == "Windows" & parallel ){
    cat("\n The parallelized version of the PIMP-algorithm are not available on Windows !! \n")
    parallel = FALSE
  }
  ##################################################################
  # parallel ?
  if(parallel) {
    # The number of cores to use, i.e. at most how many child processes will
    # be run simultaneously.
    d_ncores = parallel::detectCores()
    if(ncores>d_ncores)  stop("ncores: The number of cores is too large")

    if(ncores==0){
      ncores = max(c(1,floor(d_ncores/2)))
    }

    ## To make this reproducible
    # Random Number Generation
    # A "combined multiple-recursive generator" from L'Ecuyer (1999)
    if("L'Ecuyer-CMRG" != RNGkind()[1]) {
      cat("\n The random number generator was set to L'Ecuyer-CMRG !! \n")
      RNGkind("L'Ecuyer-CMRG")
    }
  } else { ncores = 1 }
  # specify seeds
  set.seed(seed)

  ##################################################################
  # classification ?
  classRF = is.factor(y)
  if (classRF) {
      if(rForest$type == "regression") stop("rForest$type = regression !! y a factor ")
      ##################################################################
      # permutes the response vector
      y.s=replicate(S,sample(as.integer(y)))
      # Number of permutation
      i=1:S
      varip = parallel::mclapply(i, function (i) {
                                        randomForest::randomForest(X,as.factor(y.s[,i]),mtry = mtry,
                                                                   importance=TRUE,ntree=ntree,...)[[9]][,3]
                                    },
                     mc.cores =ncores)

       varip=simplify2array(varip)
  } else {
      if(rForest$type == "classification") stop("rForest$type = classification !! y not a factor ")
      ##################################################################
      # permutes the response vector
      y.s=replicate(S,sample(y ))

      # Number of permutation
      i=1:S
      varip = parallel::mclapply(i, function (i) {
                                      randomForest::randomForest(X,y.s[,i],mtry = mtry,
                                                                 importance=TRUE,ntree=ntree,...)[[7]][,1]
                                    },
                       mc.cores =ncores)

      varip=simplify2array(varip)

  }
  #
  dimNames=dimnames(rForest$importance)[[1]]
  out=list(VarImp = matrix(rForest$importance[,1], ncol=1, dimnames=list(dimNames ,"VarImp")),
           PerVarImp = varip,
           type =   if (classRF) {"classification"} else{"regression"},
           call = match.call())


  class(out) = "PIMP"
  return(out)

  }

print.PIMP = function(x, ...){
  if (!inherits(x, "PIMP")) stop(" is not of class PIMP")
  cat("Call:\n \n")
  print(x$call)
  cat("type: ")
  print(x$type)
  cat("\noriginal VarImp:\n")
  print(x$VarImp)
  cat("\npermutation VarImp:\n")
  print(x$PerVarImp)
}
