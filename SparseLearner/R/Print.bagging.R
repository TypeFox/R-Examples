#' Print a bagging object.
#' 
#' This function prints a summary of the bagging object fitted by the bagging.lasso function.
#' 
#' @param x a fitted bagging object.
#' @export
#' @import glmnet
#' @import RankAggreg
#' @import mlbench
#' @references
#' [1] Guo, P., Zeng, F., Hu, X., Zhang, D., Zhu, S., Deng, Y., & Hao, Y. (2015). Improved Variable 
#' Selection Algorithm Using a LASSO-Type Penalty, with an Application to Assessing Hepatitis B 
#' Infection Relevant Factors in Community Residents. PLoS One, 27;10(7):e0134151.
#' [2] Tibshirani, R. (1996). Regression Shrinkage and Selection via the Lasso. Journal of the royal 
#' statistical society series B (statistical methodology), 73(3):273-282. 
#' [3] Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5-32. 
#' @examples
#' library(mlbench)
#' set.seed(0123)
#' mydata <- mlbench.threenorm(100, d=10)
#' x <- mydata$x
#' y <- mydata$classes
#' mydata <- as.data.frame(cbind(x, y))
#' colnames(mydata) <- c(paste("A", 1:10, sep=""), "y")
#' mydata$y <- ifelse(mydata$y==1, 0, 1)
#' # Split into training and testing data.
#' S1 <- as.vector(which(mydata$y==0))
#' S2 <- as.vector(which(mydata$y==1))
#' S3 <- sample(S1, ceiling(length(S1)*0.8), replace=FALSE)
#' S4 <- sample(S2, ceiling(length(S2)*0.8), replace=FALSE)
#' TrainInd <- c(S3, S4)
#' TestInd <- setdiff(1:length(mydata$y), TrainInd)
#' TrainXY <- mydata[TrainInd, ]
#' TestXY <- mydata[TestInd, ]
#' # Fit a bagging LASSO linear regression model, where the parameters 
#' # of M in the following example is set as small values to reduce the 
#' # running time, however the default value is proposed.
#' Bagging.fit <- Bagging.lasso(x=TrainXY[, -10], y=TrainXY[, 10], 
#' family=c("gaussian"), M=2, predictor.subset=round((9/10)*ncol(x)), 
#' predictor.importance=FALSE, trimmed=FALSE, weighted=TRUE, seed=0123)
#' # Print a 'bagging' object fitted by the Bagging.fit function. 
#' Print.bagging(Bagging.fit)
Print.bagging=function(x){
  if(class(x)!="bagging"){
    stop("X should be of a 'bagging' object!")
                         }
   print(x)
}