#'@title Blindfold crossvalidation of predictions from matrixpls results
#'
#'@description
#'\code{matrixpls.blindfold} Calculates blindfold predictions using \code{matrixpls}.
#'
#'@details In blindfolding, the data are first divided into n equal sized groups. Then a 
#'statistic or a model is estimated using the data after omitting one of the groups.
#'The results are then used to predict the observations in that group. The process is repeated when
#'all groups have been predicted.
#'
#'@param data Matrix or data frame containing the raw data.
#'
#'@param nGroup The number of groups to divide the data into. Setting \code{nGroup} to the number
#'of observations produces jackknife predictions.
#'
#'@param predictFun The function used to calculate the predictions.
#'
#'@param ... All other arguments are passed through to \code{\link{matrixpls}}.
#'
#'@return A matrix of class \code{matrixpls.blindfold} containing predictions calculated with 
#'blindfolding.
#'
#'@export
#'

matrixpls.blindfold <- function(data, ..., predictFun = stats::predict, nGroup = 4){
  
  
  data <- as.matrix(data)
  
  assertive::assert_all_are_in_closed_range(nGroup,1,nrow(data))
  
  groups <- ceiling(1:nrow(data)/nrow(data)*nGroup)

  
  blindfold.out <- lapply(1:nGroup,function(group){
    S <- stats::cov(data[groups!=group,])
    matrixpls.res <- matrixpls(S,...)
    predictFun(matrixpls.res,data[groups==group,])
  })
  
  
  blindfold.out <- do.call(rbind, blindfold.out)
  
  class(blindfold.out) <- c("matrixpls.blindfold", "matrix")
  attr(blindfold.out, "nGroup") <- nGroup
  blindfold.out
}

#'@title Q2	predictive relevance statistics
#'
#'@description Calculates Q2 predictive relevance statistics based on matrixpls results and 
#'crossvalidated predictions.
#'
#'@details The Q2 statistic is calculated as \code{1-sse/sso} where \code{sse} is the sum of 
#'squared prediction errors based on comparison of the \code{originalData} and 
#'\code{predictedData} and \code{sso} is based on prediction with mean.
#'
#'@param object matrixpls estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param originalData A matrix or a data.frame containing the original data.
#'
#'@param predictedData A matrix or a data.frame containing the predicted data that are compared
#'against the original data to calculate the predictive relevance statistic.
#'
#'@return A list with \code{total}, \code{block}, and \code{indicator} elements containing
#'the Q2 predictive relevance statistics for the full dataset, for each indicator block, and
#'for each indicator
#'
#'@export
#'

q2 <- function(object, originalData, predictedData){
  
  nGroup <- attr(predictedData,"nGroup")
  
  if(is.null(nGroup)) nGroup <- nrow(predictedData)
  
  # Choose only variables that exists in both data frames / matrices and 
  # sort the variables into same order
  
  v <- intersect(colnames(originalData), colnames(predictedData))
  originalData <- originalData[,v]
  predictedData <- predictedData[,v]
  
  groups <- ceiling(1:nrow(predictedData)/nrow(predictedData)*nGroup)
  
  blindfold.mean <- lapply(1:nGroup,function(group){
    means <- apply(originalData[groups!=group,v],2,mean)
    matrix(means, sum(groups!=group), length(v), byrow = TRUE)
  })
  
  blindfold.mean <- do.call(rbind, blindfold.mean)
  
  sse <- apply((originalData-predictedData)^2,2,sum)
  sso <- apply((originalData-blindfold.mean)^2,2,sum)
  
  reflective <- attr(object,"model")$reflective
  
  q2 <- list(total = 1-sum(sse)/sum(sso),
             block = apply(reflective[v,],2,function(inBlock){
               1-sum(sse[inBlock==1])/sum(sso[inBlock==1])
             }),
             indicator = 1-sse/sso)
  
  class(q2) <- "matrixplsq2"
  q2
}

#'@S3method print matrixplsq2

print.matrixplsq2 <- function(x, ...){
  cat("\n Q2 predictive relevance statistics\n")
  cat("\n Overall Q2\n")
  cat(x$total, ...)
  cat("\n Block Q2\n")
  print.table(x$block, ...)
  cat("\n Indicator Q2\n")
  print.table(x$indicator, ...)
}
