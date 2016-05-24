kfolds2Pressind <- function(pls_kfolds) {
    if (length(pls_kfolds$results_kfolds)==1) {pressind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      pressind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          pressind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }
    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                pressind_kfolds[[nnkk]][[ii]] <- (pls_kfolds$results_kfolds[[nnkk]][[ii]]-pls_kfolds$dataY_kfolds[[nnkk]][[ii]])^2
                } else {
                pressind_kfolds[[nnkk]][[ii]] <- attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$results_kfolds[[nnkk]][[ii]]-pls_kfolds$dataY_kfolds[[nnkk]][[ii]])^2
                }
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                pressind_kfolds[[nnkk]][[ii]] <- colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2)
                } else {
                pressind_kfolds[[nnkk]][[ii]] <- colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2)
                }
            }
        }
    }
rm(ii)
rm(nnkk)
return(pressind_kfolds)
}
