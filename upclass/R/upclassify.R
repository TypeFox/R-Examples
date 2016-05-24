.packageName <- 'upclass'

upclassify <-
  function (Xtrain, cltrain, Xtest, cltest = NULL, modelscope = NULL, 
            tol = 10^-5, iterlim = 1000, Aitken = TRUE, 
            ...) 
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
      res[[modelName]] <- upclassifymodel(Xtrain, cltrain, Xtest, 
                                          cltest, modelName, tol = tol, iterlim = iterlim, 
                                          Aitken = Aitken, ...)
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
