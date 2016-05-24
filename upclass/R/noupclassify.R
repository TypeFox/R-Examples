.packageName <- 'upclass'

noupclassify <-
  function (Xtrain, cltrain, Xtest, cltest = NULL, modelscope = NULL, ...) 
  {
    if (is.null(modelscope)) {
      if (is.matrix(Xtrain)) {
        d <- ncol(Xtrain)
      }
      else {
        d <- 1
      }
      modelscope <- modelvec(d)
    }
    res <- list()
    bestBIC <- -Inf
    res[["Best"]] <- list()
    for (modelName in modelscope) {
      res[[modelName]] <- list()
      res[[modelName]] <- noupclassifymodel(Xtrain, cltrain, Xtest, 
                                            cltest, modelName, ...)
      if(!is.na(res[[modelName]]$bic)){
        if (res[[modelName]]$bic > bestBIC) {
          res[["Best"]] <- res[[modelName]]
          bestBIC <- res[[modelName]]$bic
        }
      }
    }
    class(res)<-"upclassfit" 
    res
  }
