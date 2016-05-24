compVarImp = function(X, y,rForest,nPerm = 1){

  if (!inherits(rForest, "randomForest")) stop("rForest is not of class randomForest")
  if(is.null(rForest$inbag)) stop("in randomForest keep.inbag must be set to True")
  if(is.null(rForest$forest)) stop("in randomForest keep.forest must be set to True")

  n = nrow(X)
  p = ncol(X)
  if (n == 0) stop("data (x) has 0 rows")
  x.col.names <- if (is.null(colnames(X))) 1:p else colnames(X)

  if (is.data.frame(X)) {
      X = data.matrix(X)
  }

  if (!is.null(y)) {
    if (length(y) != n) stop("length of response must be the same as predictors")
  }

  if (any(is.na(X))) stop("NA not permitted in predictors")
  if (any(is.na(y))) stop("NA not permitted in response")

  classRF = is.factor(y)
  if (classRF) {
      if (!all(levels(y) == levels(rForest$y)))  stop("y and rForest$y must have the same levels")
      if(rForest$type == "regression") stop("rForest$type = regression !! y a factor ")
      nclass = length(levels(y))
      classOut = .Call('vita_Rcpp_compVarImpCL', PACKAGE = 'vita',
                       t(X), as.integer(y), n, p,
                       rForest$ntree, as.integer(nclass),
                       rForest$forest$treemap[,1,], rForest$forest$treemap[,2,],
                       rForest$forest$nodestatus, rForest$forest$xbestsplit,
                       rForest$forest$nodepred, rForest$forest$bestvar,
                       rForest$inbag, rForest$forest$ndbigtree,
                       rForest$forest$ncat, rForest$forest$maxcat )

      dimnames(classOut$importance)   = list(x.col.names, c(levels(y), "MeanDecreaseAccuracy"))
      dimnames(classOut$importanceSD) = list(x.col.names, c(levels(y), "MeanDecreaseAccuracy"))
      return( list( importance = classOut$importance,
                    importanceSD = classOut$importanceSD,
                    type = "classification"
                  )
              )

  }else {
      if(rForest$type == "classification") stop("rForest$type = classification !! y not a factor ")

      regOut = .Call('vita_Rcpp_compVarImpReg', PACKAGE = 'vita', t(X),y, n,p,
                   rForest$ntree,nPerm,rForest$forest$leftDaughter,
                   rForest$forest$rightDaughter,rForest$forest$nodestatus,
                   rForest$forest$xbestsplit,rForest$forest$nodepred,
                   rForest$forest$bestvar,rForest$inbag,rForest$forest$ndbigtree,
                   rForest$forest$ncat, max(rForest$forest$ncat))

      return( list( importance =  matrix(regOut$importance, p, 1,
                                         dimnames=list(x.col.names,  c("%IncMSE")) ),
                    importanceSD = matrix(regOut$importanceSD, p, 1,
                                          dimnames=list(x.col.names,  c("%IncMSE")) ),
                    type = "regression"
                  )
              )

  }


}
