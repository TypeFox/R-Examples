VarImpCVl = function (X_l, y_l, rForest, nPerm = 1){

  if (!inherits(rForest, "randomForest")) stop("rForest is not of class randomForest")
  if(is.null(rForest$forest)) stop("in randomForest keep.forest must be set to True")

  n = nrow(X_l)
  p = ncol(X_l)
  if(n == 0) stop("data (X_l) has 0 rows")
  x.col.names <- if (is.null(colnames(X_l))) 1:p else colnames(X_l)

  if (is.data.frame(X_l)) {
    xlevels = lapply(X_l, function(x) if (is.factor(x)) levels(x) else 0)
    ncat = sapply(xlevels, length)
    ## Treat ordered factors as numerics.
    ncat = ifelse(sapply(X_l, is.ordered), 1, ncat)
    maxncat = max(ncat)
    X_l = data.matrix(X_l)

  } else{
    ncat = rep(1, p)
    maxncat = 1
  }

  if (!is.null(y_l)) {
    if (length(y_l) != n) stop("length of response must be the same as predictors")
  }

  if (any(is.na(X_l))) stop("NA not permitted in predictors")
  if (any(is.na(y_l))) stop("NA not permitted in response")

  classRF = is.factor(y_l)
  # classification ?
  if (classRF) {

    if (!all(levels(y_l) == levels(rForest$y_l)))  stop("y_l and rForest$y must have the same levels")
    if(rForest$type == "regression") stop("rForest$type = regression !! y_l a factor ")
    nclass = length(levels(y_l))

    classOut = .Call('vita_Rcpp_VarImpCVLCL', PACKAGE = 'vita',
                     t(X_l), as.integer(y_l), n, p,
                     rForest$ntree, as.integer(nclass),
                     rForest$forest$treemap[,1,], rForest$forest$treemap[,2,],
                     rForest$forest$nodestatus, rForest$forest$xbestsplit,
                     rForest$forest$nodepred, rForest$forest$bestvar,
                     rForest$forest$ndbigtree, ncat, maxncat)

    return(  list( fold_importance  =  matrix(classOut$fold_importance, p, 1,
                                               dimnames=list(x.col.names,  c("MeanDecreaseAccuracy")) ),
                   type = "classification"
                   )
             )

    } else{
      if(rForest$type == "classification") stop("rForest$type = classification !! y_l not a factor ")

      regOut = .Call('vita_Rcpp_VarImpCVLReg', PACKAGE = 'vita', t(X_l),y_l, n, p,
                     rForest$ntree,nPerm, rForest$forest$leftDaughter,
                     rForest$forest$rightDaughter, rForest$forest$nodestatus,
                     rForest$forest$xbestsplit, rForest$forest$nodepred,
                     rForest$forest$bestvar, rForest$forest$ndbigtree, ncat, maxncat)

      return( list( fold_importance =  matrix(regOut$fold_importance, p, 1,
                                         dimnames=list(x.col.names,  c("%IncMSE")) ),
                    type = "regression"
                   )
              )
    }

}
