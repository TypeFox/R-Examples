#' @family performance analyzers
#' @export
#' 
#' @title
#' Perform generic out-of-bag error analysis.
#' @description
#' If performing regression, calculate which out-of-bag residuals and MSE.
#' Otherwise, calculate which out-of-bag observations were classified correctly,
#' what the overall misclassification rate is, as well as the confusion matrix. 
#' 
#' @param prediction a vector of predicted responses.
#' @param response a vector of true response.
#' @param oobObs a vector of indices which values in \code{predictions} are of
#' out-of-bag observations. 
#' 
#' @return
#' If performing regression, return a list with components:
#' \item{oobMSE}{the out-of-bag mean squared error.}
#' \item{resVec}{a vector of length \code{nrow(data)} whose entries correspond to
#' observations in \code{data}. The entry has values \code{NA} if the observation
#' was not out-of-bag, and the difference between the predicted and true response
#' (the residual) if the observation was out-of-bag.}
#' 
#' Otherwise, return a list with components:
#' \item{oobErr}{overall misclassification rate.}
#' \item{oobConfMat}{the confusion matrix of out-of-bag predictions against the 
#' true class labels.}
#' \item{errVec}{a vector of length \code{nrow(data)} whose entries correspond to
#' observations in \code{data}. The entry has values \code{NA} if the observation
#' was not out-of-bag, and a 1 or 0 depending whether \code{estimator} failed to
#' correctly classify the observation.}


defaultOOBPerformanceAnalysis <- function(prediction, response, oobObs) {
  
  n <- nrow(prediction)
  oobPreds <- prediction[oobObs]
  oobResponse <- response[oobObs]
  
  if (class(oobResponse) %in% c("factor", "character")) {
    oobPreds <- as.character(oobPreds)
    oobResponse <- as.character(oobResponse)
    
    errVec <- rep.int(NA, length(prediction))
    errVec[oobObs] <- as.numeric(oobPreds != oobResponse)
    
    oobConfMat <- table(oobPreds, oobResponse)
    oobErr <- mean(errVec, na.rm=TRUE)
    
    list(oobErr=oobErr, oobConfMat=oobConfMat, errVec=errVec)
    
  } else {

    resVec <- rep.int(NA, n)
    resVec[oobObs] <- oobPreds - oobResponse
    
    oobMSE <- mean(resVec^2, na.rm=TRUE)
    
    list(oobMSE=oobMSE, resVec=resVec)
  }
}

class(defaultOOBPerformanceAnalysis) <- c("performanceAnalyzer",
                                          class(defaultOOBPerformanceAnalysis))