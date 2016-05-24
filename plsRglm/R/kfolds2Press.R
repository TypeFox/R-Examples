kfolds2Press <- function(pls_kfolds) {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
      max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],ncol)))
      press_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      press_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],ncol)))
          press_kfolds[[jj]] <- rep(0,max_nt[jj])
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
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]]-pls_kfolds$dataY_kfolds[[nnkk]][[ii]])^2
                } else {
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]]-pls_kfolds$dataY_kfolds[[nnkk]][[ii]])^2
                }
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2) 
                } else {
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2)           
                }
            }
        }
    }
rm(ii)
rm(nnkk)
return(press_kfolds)
}
