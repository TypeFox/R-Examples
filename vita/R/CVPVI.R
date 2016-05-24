CVPVI = function(X, ...)  UseMethod("CVPVI")

CVPVI.default = function(X, y, k = 2,
                         mtry= if (!is.null(y) && !is.factor(y))
                               max(floor(ncol(X)/3), 1) else floor(sqrt(ncol(X))),
                         ntree = 500, nPerm = 1, parallel = FALSE,
                         ncores = 0, seed = 123, ...){

  ##################################################################
  n = nrow(X)
  p = ncol(X)
  if (n == 0) stop("data (x) has 0 rows")
  x.col.names <- if (is.null(colnames(X))) 1:p else colnames(X)


  if (!is.null(y)) {
    if (length(y) != n) stop("length of response must be the same as predictors")
  }

  if (any(is.na(X))) stop("NA not permitted in predictors")
  if (any(is.na(y))) stop("NA not permitted in response")
  ##################################################################
  # Windows and parallel ?
  if(Sys.info()[['sysname']] == "Windows" & parallel ){
    cat("\n The parallelized version of the CVPVI implementation are not available on Windows !! \n")
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

      if(k < ncores ) {
        ncores2 = floor(ncores/k)
        ncores = k
      } else {
        ncores2 = 1
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
  ## Split indexes k- folds
  cuts = round(length(y)/k)
  from = (0:(k-1)*cuts)+1
  to = (1:k*cuts)
  if(to[k]!=length(y)) to[k] = length(y)
  rs = sample(1:length(y))
  l = 1:k

  ##################################################################
  #  parallel CVPVI implementation
  if(parallel) {

      ##################################################################
      ## Split ntree/ncores2
      nt = rep(floor(ntree/ncores2),ncores2)
      if(sum(nt)<ntree){nt[ncores2] = nt[ncores2] + ntree - sum(nt)}

      cvl_varim = parallel::mclapply(l, function (l) {
                                    lth = rs[from[l]:to[l]]
                                    # without the l-th data set
                                    Xl = X[-lth,]
                                    yl = y[-lth]
                                    # the l-th data set
                                    X_l = X[lth,]
                                    y_l = y[lth]

                                    CVVI_ln = parallel::mclapply(nt, function (nt) {
                                                     rForest_ln = randomForest::randomForest(Xl,yl,
                                                                                       mtry = mtry,
                                                                                       ntree = nt,
                                                                                       keep.forest = T,...)

                                                     VarImpCVl(X_l,y_l,rForest_ln)[[1]]
                                                },mc.cores = ncores2)

                                    CVVI_ln=simplify2array(CVVI_ln)[,1,]
                                    return((CVVI_ln%*%nt)/ntree)

                   }, mc.cores = ncores)

  ##################################################################
  # non parallel CVPVI implementation
  } else {
          nt=1
          cvl_varim = lapply(l, function (l) {
                                    lth = rs[from[l]:to[l]]
                                    # without the l-th data set
                                    Xl = X[-lth,]
                                    yl = y[-lth]
                                    # the l-th data set
                                    X_l = X[lth,]
                                    y_l = y[lth]

                                    rForest_l = randomForest::randomForest(Xl,yl,
                                                                           mtry = mtry,
                                                                           ntree = ntree,
                                                                           keep.forest = T,...)
                                    # Compute l-th fold-specific variable importance
                                    return(VarImpCVl(X_l,y_l,rForest_l,nPerm))
                              }
                      )

  }

  if (length(nt) == 1 ) {
   m = as.data.frame.list(cvl_varim)
   cvl_varimf = m[seq(from = 1, to = 2*k, by = 2)]
   rownames(cvl_varimf) = x.col.names
   colnames(cvl_varimf) = paste0(l, rep("-fold_PerVarImp", k))
  } else {
    cvl_varimf = simplify2array(cvl_varim)[,1,]
    rownames(cvl_varimf) = x.col.names
    colnames(cvl_varimf) = paste0(l, rep("-fold_PerVarImp", k))
  }
 out = list(
          # fold-specific permutation variable importance
          fold_varim = cvl_varimf,
          # cross-valdated permutation variable importance
          cv_varim = matrix(rowMeans(cvl_varimf),ncol = 1,dimnames = list(x.col.names ,"CV_PerVarImp")),
          type = if (is.factor(y)) {"classification"} else{"regression"},
          call = match.call())
 class(out) = "CVPVI"
 return(out)

}


print.CVPVI = function(x, ...){
  if (!inherits(x, "CVPVI")) stop(" is not of class CVPVI")
  cat("Call:\n \n")
  print(x$call)
  cat("type: ")
  print(x$type)
  cat("\nfold-specific permutation variable importance:\n\n")
  print(x$fold_varim)
  cat("\ncross-validated permutation variable importance:\n\n")
  print(x$cv_varim)
}
