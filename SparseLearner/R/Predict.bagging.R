#' Make predictions for new data from a 'bagging' object.
#' 
#' This function makes predictions for new data from a bagging LASSO linear or logistic regression model, using the stored 'bagging' object, with or without the use of trimmed bagging strategy. 
#' 
#' @param object a fitted 'bagging' object.
#' @param newx matrix of new values for x at which predictions are to be made. Must be a matrix. See documentation for Bagging.lasso.
#' @param y response variable. Defaults to NULL. If the response variable for the newx matrix is known and input, the corresponding validation measures can be calculated for evaluating prediction performance.  
#' @param trimmed logical. Should a trimmed bagging strategy be performed? Defaults to FALSE. This argument should correspond to the same setting in the Bagging.lasso function. See documentation for Bagging.lasso.
#' @param scale.trimmed the portion to trim of the "worst" based-level models, in the sense of having the largest error rates, and to average only over the most accurate base-level models. Defaults to 0.75. 
#' @export
#' @import glmnet
#' @import RankAggreg
#' @import mlbench
#' @references
#' [1] Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5-32.
#' [2] Croux, C., Joossens, K., & Lemmens, A. (2007). Trimmed bagging. Computational Statistics & 
#' Data Analysis, 52(1), 362-368. 
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
#' Bagging.fit
#' # Make predictions from a bagging LASSO linear regression model.
#' pred <- Predict.bagging(Bagging.fit, newx=TestXY[, -10], y=NULL)
#' pred
Predict.bagging <- function(object, newx, y=NULL, trimmed=FALSE, scale.trimmed=0.75){
rmse <- function(truth, predicted){
    predicted <- as.numeric(predicted)
    mse <- mean((truth-predicted)*(truth-predicted))
    rmse  <- sqrt(mse)
    rmse
                                  }
mae <- function(truth, predicted){
    predicted <- as.numeric(predicted)
    mae <- mean(abs(predicted-truth))
    mae
                                 }
re <- function(truth, predicted){
    predicted <- as.numeric(predicted)
    mse <- mean(((truth-predicted)/truth)*((truth-predicted)/truth))
    re  <- sqrt(mse)
    re
                                }
smape <- function(truth, predicted){
    predicted <- as.numeric(predicted)
    smape <- mean(abs(truth-predicted)/((abs(truth)+abs(predicted))/2))
    smape
                                   }
accuracy <- function(truth, predicted){
    if(length(truth) > 0)
        sum(truth==predicted)/length(truth) else    
        return(0)
                                      }
sensitivity <- function(truth, predicted){
    # 1 means positive (present)
    if(sum(truth==1) > 0)
        sum(predicted[truth==1]==1)/sum(truth==1) else
        return(0)
                                          }
specificity <- function(truth, predicted){
    if(sum(truth==0) > 0)
        sum(predicted[truth==0]==0)/sum(truth==0) else
        return(0)
                                         }
AUC <- function(truth, probs){
    # probs - probability of class 1
    q <- seq(0, 1, .01)
    sens <- rep(0, length(q))
    spec <- rep(0, length(q))
    ly <- levels(truth)
    for(i in 1:length(q)){
        pred <- probs >= q[i]
        pred[pred] <- 1
        pred <- factor(pred, levels=ly)
        sens[i] <- sensitivity(truth, pred)
        spec[i] <- specificity(truth, pred)
                         }
    # make sure it starts and ends at 0, 1
    sens <- c(1, sens, 0)
    spec <- c(0, spec, 1)
    trap.rule <- function(x,y) sum(diff(x)*(y[-1]+y[-length(y)]))/2 
    auc <- trap.rule(rev(1-spec), rev(sens))
    auc
                             }
kia <- function(truth, predicted){
    TP   <- sum(predicted[truth==1]==1)
    TN   <- sum(predicted[truth==0]==0)
    FN   <- sum(truth==1)-TP
    FP   <- sum(truth==0)-TN
    N    <- TP+TN+FN+FP
    Pobs <- (TP+TN)/N
    Pexp <- ((TP+FN)*(TP+FP)+(FP+TN)*(FN+TN))/N^2
    kia  <- (Pobs-Pexp)/(1-Pexp)
    kia
                                 }
convertScores <- function(scores){
    scores <- t(scores)
    ranks <- matrix(0, nrow(scores), ncol(scores))
    weights <- ranks
    for(i in 1:nrow(scores)){
        ms <- sort(scores[i,], decreasing=TRUE, ind=TRUE)
        ranks[i,] <- colnames(scores)[ms$ix]
        weights[i,] <- ms$x
                            }
    list(ranks = ranks, weights = weights)
                                 }

  Model.type <- object$family
  if (Model.type==c("gaussian")){
    Models <- object$models.fitted
    M <- object$M
    n <- nrow(newx)
    predicted <- matrix(0, n, M)
    model.averaging <- c("mean")
    if(trimmed){
        if(is.null(object$models.trimmed)) { stop("Trimming is not performed!") }
        M <- ceiling(M*scale.trimmed)
        Models <- object$models.trimmed
        predicted <- matrix(0, n, M)
               }

    for(i in 1:M){
        cvfit1 <- Models[[i]]
        model.final <- cvfit1$glmnet.fit
        s2 <- rownames(as.matrix(coef(model.final, s=cvfit1$lambda.min)))
        testing <- as.matrix(newx[, s2[-1]]) 
        predicted[, i] <- predict(cvfit1, newx=testing, s=c("lambda.min"), type=c("link"))
                 }
     predicted <- matrix(as.numeric(predicted), n, M)

       if(model.averaging==c("mean")){
           newvalue <- apply(predicted, 1, mean)
           sevalue <- apply(predicted, 1, function(x) sd(x)/sqrt(length(x)))
                                     }
    res <- list()
    if(!is.null(y)){
        valM  <- c("rmse", "mae", "re", "smape")
        rmse  <- rmse(y, newvalue)
        mae   <- mae(y, newvalue)
        re    <- re(y, newvalue)
        smape <- smape(y, newvalue)
        bagging.prediction <- matrix(c(rmse, mae, re, smape), 1, 4)
        colnames(bagging.prediction) <- valM
        rownames(bagging.prediction) <- "bagging"
                    }
    if(is.null(y))
        res <- list(y.new=newvalue, y.se=sevalue, predicted.matrix=predicted)
    else
        res <- list(y.new=newvalue, y.se=sevalue, predicted.matrix=predicted, bagging.prediction=bagging.prediction)
                               } # end of linear regression model


  if (Model.type==c("binomial")){
    Models <- object$models.fitted
    M <- object$M
    n <- nrow(newx)
    predicted <- matrix(0, n, M)
    model.averaging <- c("majorvoting")
    if(trimmed){
        if(is.null(object$models.trimmed)) { stop("Trimming is not performed!") }
        M <- ceiling(M*scale.trimmed)
        Models <- object$models.trimmed
        predicted <- matrix(0, n, M)
               }

    for(i in 1:M){
        cvfit1 <- Models[[i]]
        model.final <- cvfit1$glmnet.fit
        s2 <- rownames(as.matrix(coef(model.final, s=cvfit1$lambda.min)))
        testing <- as.matrix(newx[, s2[-1]]) 
        predicted[, i] <- predict(cvfit1, newx=testing, s=c("lambda.min"), type=c("class"))
                 }
     predicted <- matrix(as.numeric(predicted), n, M)

       if(model.averaging==c("majorvoting")){
           new.class <- factor(apply(predicted, 1, function(x) ifelse(sum(x) > floor(M/2), 1, 0)), levels=c("0", "1"))
           probabilities <- apply(predicted, 1, function(x) sum(x)/M)
                                             }
    res <- list()
    if(!is.null(y)){
        valM <- c("accuracy", "sensitivity", "specificity", "auc", "kia")
        acc  <- accuracy(y, new.class)
        sens <- sensitivity(y, new.class)
        spec <- specificity(y, new.class)
        kia  <- kia(y, new.class)
        auc  <- AUC(y, probabilities)     
        bagging.prediction <- matrix(c(acc, sens, spec, auc), 1, 4)
        colnames(bagging.prediction) <- valM
        rownames(bagging.prediction) <- "bagging"
                    }
    if(is.null(y))
        res <- list(y.new=new.class, probabilities=probabilities, predicted.matrix=predicted)
    else
        res <- list(y.new=new.class, probabilities=probabilities, predicted.matrix=predicted, bagging.prediction=bagging.prediction)
                                  }# end of logistic regression model


    class(res) <- "BaggingPrediction"
    res
}
