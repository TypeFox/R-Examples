#' Generate a plot of variable importance.
#' 
#' This function generates a plot for evaluating variable importance based on a bagging object fitted by the bagging.lasso model.
#' 
#' @param x a fitted bagging object.
#' @param max.var.show the maximum number of variables to be shown in the plot. Defaults to 40.
#' @param xlab a title for the x axis. 
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot. 
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
#' predictor.importance=TRUE, trimmed=FALSE, weighted=TRUE, seed=0123)
#' Plot.importance(Bagging.fit)
Plot.importance=function(x, max.var.show=40, xlab="Importance Score", ylab=NULL, main="Variable Importance Plot"){
    if(class(x)!="bagging"){
      stop("X should be of a 'bagging' object!")
                         }
    varImp <- abs(x$importance[, 1])
    vars <- order(varImp, decreasing=TRUE)
    p <- length(varImp)
    pm <- names(varImp)
    max.v <- max.var.show
    if(p<max.v)
       max.v <- p
    ylab=ylab
    xlab=xlab
    main=main
    dotchart(varImp[vars[max.v:1]], labels=pm[vars[max.v:1]], ylab=ylab, xlab=xlab, main=main)
}
