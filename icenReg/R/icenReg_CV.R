abs_inv <- function(pred, t_val){
 mean(abs(1/(pred+1) - 1/(t_val+1 ) ) ) 
}

evalCV <- function(imputes, preds, cv_fun){
  ans <- numeric(2)
  names(ans) <- c('mean loss', 'imputation SE')
  nImputes <- ncol(imputes)
  if(nImputes == 0) {ans[0] <- NA; return(ans)}
  lossFunEsts <- numeric(nImputes)
  for(i in 1:nImputes)
    lossFunEsts[i] <- cv_fun(preds, imputes[,i])
  
  ans[1] <- mean(lossFunEsts)
  ans[2] <- sd(lossFunEsts) / sqrt(nImputes)
  return(ans)
}

icenReg_cv <- function(fit, loss_fun = abs_inv, folds = 10, numImputes = 100, useMCore = F){
  if(folds == 1) stop('folds must be greater than 1')
  
  if(useMCore) `%myDo%` <- `%dopar%`
  else `%myDo%` <- `%do%`
  
  numRow <- nrow(fit$getRawData() )
  removePerFold <- floor(numRow/folds)
  inds <- 1:numRow
  sampInds <- sample(inds, numRow)
  grabInds <- 1:removePerFold
  cv_inds <- list()
  for(i in 1:(folds-1) ){
    cv_inds[[i]] <- sampInds[1:removePerFold]
    sampInds <- sampInds[-(1:removePerFold)]
  }
  if(length(sampInds) > 0)  cv_inds[[folds]] <- sampInds
  modCall <- fit$call
  modCall$data <- as.name("TRAIN_DATA_ICENREG")
#  if(is(fit, 'sp_fit')) modCall$bs_samples = 0
  rawData <- fit$getRawData()
  cvItems <- list(data = rawData, cv_inds = cv_inds, 
                  loss_fun = loss_fun, numImputes = numImputes,
                  call = modCall)
   cvSummary <- foreach(i = 1:folds, .combine = rbind ) %myDo%
              {
              TRAIN_DATA_ICENREG <- rawData[-cv_inds[[i]], ]
              VALID_DATA_ICENREG <- rawData[cv_inds[[i]], ]
              cv_fit <- eval(modCall)
              cv_preds <- predict.icenReg_fit(cv_fit, newdata = VALID_DATA_ICENREG)
              cv_imputes <- imputeCens(cv_fit, VALID_DATA_ICENREG, numImputes = numImputes)
              ans <- evalCV(cv_imputes, cv_preds, loss_fun)
              rm(cv_fit, cv_preds, cv_imputes, TRAIN_DATA_ICENREG, VALID_DATA_ICENREG)
    return(ans)
    }
  mean_cv_error <- mean(cvSummary[,1])
  imputed_cv_se <- mean(cvSummary[,2]) / sqrt(nrow(cvSummary)) 
  ans <- c(mean_cv_error, imputed_cv_se)
  names(ans) <- c('mean_cv_error', 'imputed_cv_se')
  return(ans)
}



