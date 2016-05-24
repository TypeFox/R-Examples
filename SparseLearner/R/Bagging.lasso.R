#' A Bagging Prediction Model Using LASSO Selection Algorithm.
#' 
#' This function performs a bagging prediction for linear and logistic regression model using the LASSO selection algorithm.
#' 
#' @param x input matrix. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables.
#' @param y response variable. For family="gaussian", y is a vector of quantitative response. For family="binomial" should be a factor with two levels '0' and '1' and the level of '1' is the target class. 
#' @param family response type (see above).
#' @param M the number of base-level models (LASSO linear or logistic regression models) to obtain a final prediction. Note that it also corresponds to the number of bootstrap samples to draw. Defaults to 100.
#' @param subspace.size the number of random subspaces to construct an ensemble prediction model. Defaults to 10.
#' @param predictor.subset the subset of randomly selected predictors from the training set to reduce the original p-dimensional feature space. Defaults to (9/10)*ncol(x) where ncol(x) represents the the original p-dimensional feature space of input matrix x.
#' @param boot.scale the scale of sample size in each bootstrap re-sampling, relative to the original sample size. Defaults to 1.0, equaling to the original size of training samples.
#' @param kfold the number of folds of cross validation - default is 10. Although kfold can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is kfold=3.
#' @param predictor.importance logical. Should the importance of each predictor in the bagging LASSO model be evaluated? Defaults to TRUE. A permutation-based variable importance measure estimated by the out-of-bag error rate is adapted for the bagging model. 
#' @param trimmed logical. Should a trimmed bagging strategy be performed? Defaults to FALSE. Traditional bagging draws bootstrap samples from the training sample, applies the base-level model to each bootstrap sample, and then averages over all obtained prediction rules. The idea of trimmed bagging is to exclude the bootstrapped prediction rules that yield the highest error rates and to aggregate over the remaining ones.
#' @param weighted logical. Should a weighted rank aggregation procedure be performed? Defaults to TRUE. This procedure uses a Monte Carlo cross-entropy algorithm combining the ranks of a set of based-level model under consideration via a weighted aggregation that optimizes a distance criterion to determine the best performance base-level model. 
#' @param verbose logical. Should the iterative process information of bagging model be presented? Defaults to TRUE. 
#' @param seed the seed for random sampling, with the default value 0123.
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
#' # Example 1: Bagging LASSO linear regression model.
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
#' # Print a 'bagging' object fitted by the Bagging.fit function. 
#' Print.bagging(Bagging.fit)
#' # Make predictions from a bagging LASSO linear regression model.
#' pred <- Predict.bagging(Bagging.fit, newx=TestXY[, -10], y=NULL, trimmed=FALSE)
#' pred
#' # Generate the plot of variable importance. 
#' Plot.importance(Bagging.fit)

#' # Example 2: Bagging LASSO logistic regression model.
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
#' # Fit a bagging LASSO logistic regression model, where the parameters 
#' # of M in the following example is set as small values to reduce the 
#' # running time, however the default value is proposed.
#' Bagging.fit <- Bagging.lasso(x=TrainXY[, -11], y=TrainXY[, 11], 
#' family=c("binomial"), M=2, predictor.subset=round((9/10)*ncol(x)), 
#' predictor.importance=TRUE, trimmed=FALSE, weighted=TRUE, seed=0123)
#' # Print a 'bagging' object fitted by the Bagging.fit function. 
#' Print.bagging(Bagging.fit)
#' # Make predictions from a bagging LASSO logistic regression model.
#' pred <- Predict.bagging(Bagging.fit, newx=TestXY[, -11], y=NULL, trimmed=FALSE)
#' pred
#' # Generate the plot of variable importance. 
#' Plot.importance(Bagging.fit)
Bagging.lasso <- function(x, y, family=c("gaussian", "binomial"), M=100, subspace.size=10, predictor.subset=round((9/10)*ncol(x)), boot.scale=1.0, kfold=10, 
                          predictor.importance=TRUE, trimmed=FALSE, weighted=TRUE, verbose=TRUE, seed=0123){

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

 if (family==c("gaussian")) { 
    x <- as.matrix(x)
    y <- as.numeric(y)
    if(!is.null(seed)) {
        set.seed(seed) } else
        set.seed(0123)
    validation <- c("rmse", "mae", "re", "smape")
    distance <- c("Spearman")
    rownames(x) <- NULL
    n <- length(y)
    nvm <- length(validation)
    fittedModels <- list()
    trimmedModels <- list()
    x_length <- ncol(x)
    RecordM <- matrix(0, M, ncol(x))
    colnames(RecordM) <- colnames(x)
    model.rmse <- c()
    for(k in 1:M){
        s <- sample(round(boot.scale*n), replace=TRUE)  # Size% of original samples
        training <- x[s, ]
        testing  <- x[-unique(s), ]
        trainY <- y[s]
        # Random subspace
        Res <- list()  
        predicted <- list()
        for(n0 in 1:subspace.size){
            s0 <- sample(x=dim(x)[2], size=predictor.subset, replace=FALSE)
            training_1 <- training[,s0]
            testing_1 <- testing[,s0]
            Res[[n0]] <- cv.glmnet(x=as.matrix(training_1), y=trainY, type.measure="mse", nfolds=kfold, family="gaussian")
            predicted[[n0]] <- predict(Res[[n0]], newx=as.matrix(testing_1), s=c("lambda.min"), type=c("link"))
                                  }
        # Compute validation measures
        scores <- matrix(0, subspace.size, nvm)
        rownames(scores) <- 1:subspace.size
        colnames(scores) <- validation
        truth <- y[-unique(s)]
        for(i in 1:subspace.size){
            for(j in 1:nvm){
                scores[i,j] <- switch(validation[j],
                    "rmse"    = 1/rmse(truth, predicted[[i]]),
                    "mae"     = 1/mae(truth, predicted[[i]]),
                    "re"      = 1/re(truth, predicted[[i]]),
                    "smape"   = 1/smape(truth, predicted[[i]])
                                      )
                           }
                           }
        # Perform rank aggregation
        algorithms <- as.character(1:subspace.size)
        convScores <- convertScores(scores)
        if(nvm > 1 & subspace.size <= 8)
            if(weighted)
                fittedModels[[k]] <- Res[[which(algorithms == BruteAggreg(convScores$ranks,
                                     subspace.size, convScores$weights, distance=distance)$top.list[1])]]
            else
                fittedModels[[k]] <- Res[[which(algorithms == BruteAggreg(convScores$ranks, subspace.size,
                                     distance=distance)$top.list[1])]]
        else if(nvm > 1 & subspace.size > 8)
            if(weighted)
                fittedModels[[k]] <- Res[[which(algorithms == RankAggreg(convScores$ranks,
                                     subspace.size, convScores$weights, distance=distance, verbose=FALSE)$top.list[1])]]
            else
                fittedModels[[k]] <- Res[[which(algorithms == RankAggreg(convScores$ranks, subspace.size,
                                     distance=distance, verbose=FALSE)$top.list[1])]]        
        else
            fittedModels[[k]] <- Res[[which.max(scores[,1])]]
	 # Variable importance evaluation
		if(predictor.importance){
                   model.final <- fittedModels[[k]]$glmnet.fit
                   LassoM.coef <- coef(model.final, s=fittedModels[[k]]$lambda.min)
                   Var_subset <- names(LassoM.coef[as.vector(LassoM.coef[,1]!=0),])
                   Var_subset <- Var_subset[Var_subset!=c("(Intercept)")]

                   if(!is.null(Var_subset)){
                     training_2 <- training[, Var_subset]
                     Res1 <- cv.glmnet(x=as.matrix(training_2), y=trainY, type.measure="mse", nfolds=kfold, family="gaussian")
                     model.final1 <- Res1$glmnet.fit
                     LassoM.coef1 <- coef(model.final1, s=Res1$lambda.min)[-1,]
                        if (length(LassoM.coef1)!=0){
                          for(i0 in 1:length(LassoM.coef1)){
                            for(j0 in 1:x_length){
                              if (names(LassoM.coef1)[i0]==colnames(RecordM)[j0]){RecordM[k,j0]=LassoM.coef1[i0]}
                              else {RecordM[k,j0]=RecordM[k,j0]}
                                                  }
                                                            }
                                                     }
                                           }
		          }
        # Trimmed bagging
		if(trimmed){
                glmFit <- fittedModels[[k]]
                model.glmFit <- glmFit$glmnet.fit
                model.var <- rownames(as.matrix(coef(model.glmFit, s=glmFit$lambda.min)))
                model.testing <- as.matrix(testing[, model.var[-1]]) 
                model.predicted <- predict(glmFit, newx=model.testing, s=c("lambda.min"), type=c("link"))
                model.rmse[k] <- rmse(truth, model.predicted)
		           }
	 # Running message
        if(verbose)
            cat("Iter ", k, "\n")

                 } # Loop end 1:M

    # Variable importance socres
    RecordM <- abs(RecordM)
    varImportance <- abs(as.matrix(apply(RecordM, 2, mean), ncol(RecordM), 1))
    # Trimmed bagging
    if(trimmed){
        trimmedModels <- fittedModels
          for(r in 1:M){
            trimmedModels[[rank(model.rmse)[r]]] <- fittedModels[[r]]
		           }
               }

                           } # linear regression model end

 if (family==c("binomial")) { 
    x <- as.matrix(x)
    y <- as.factor(y)    
    ly <- levels(y)
    if(length(ly) == 2 && any(ly != c("0", "1"))){
        stop("For logistic regression model, levels in y must be 0 and 1")}
    if(!is.null(seed)) {
        set.seed(seed)} else
        set.seed(0123)
    validation <- c("accuracy", "sensitivity", "specificity", "auc", "kia")
    distance <- c("Spearman")
    rownames(x) <- NULL
    n <- length(y)
    nvm <- length(validation)
    fittedModels <- list()
    trimmedModels <- list()
    RecordM <- matrix(0, M, ncol(x))  
    colnames(RecordM) <- colnames(x)
    model.accuracy <- c()
    for(k in 1:M){
        repeat{ 
            s <- sample(round(boot.scale*n), replace=TRUE) # Size% of original samples
            if(length(table(y[s])) >= 2 & length(table(y[-s])) >= 2)
                break
              }
        training <- x[s, ]
        testing  <- x[-unique(s), ]
        trainY <- y[s]
        # Random subspace
        Res <- list()  
        probabilities <- list()
        predicted <- list()
        for(n0 in 1:subspace.size){
            s0 <- sample(length(colnames(x)), replace=FALSE)
            s0 <- s0[1:predictor.subset]
            training_1 <- training[, s0]
            testing_1 <- testing[, s0]
            Res[[n0]] <- cv.glmnet(x=as.matrix(training_1), y=as.factor(trainY), type.measure="deviance", nfolds=kfold, family="binomial")
            predicted[[n0]] <- predict(Res[[n0]], newx=as.matrix(testing_1), s=c("lambda.min"), type=c("class"))
            probabilities[[n0]] <- predict(Res[[n0]], newx=as.matrix(testing_1), s=c("lambda.min"), type=c("response"))
                                  }
        # Compute validation measures           
        scores <- matrix(0, subspace.size, nvm)
        rownames(scores) <- 1:subspace.size
        colnames(scores) <- validation
        truth <- y[-unique(s)]
        for(i in 1:subspace.size){
            for(j in 1:nvm){
                scores[i,j] <- switch(validation[j],
                    "accuracy"    = accuracy(truth, factor(predicted[[i]], levels=ly)),
                    "sensitivity" = sensitivity(truth, factor(predicted[[i]], levels=ly)),
                    "specificity" = specificity(truth, factor(predicted[[i]], levels=ly)),
                    "kia"         = kia(truth, factor(predicted[[i]], levels=ly)),
                    "auc"         = AUC(truth, probabilities[[i]])
                                     )     
                           }
                           }
        # Perform rank aggregation
        algorithms <- as.character(1:subspace.size)
        convScores <- convertScores(scores)
        if(nvm > 1 & subspace.size <= 8)
            if(weighted)
                fittedModels[[k]] <- Res[[which(algorithms == BruteAggreg(convScores$ranks,
                                     subspace.size, convScores$weights, distance=distance)$top.list[1])]]
            else
                fittedModels[[k]] <- Res[[which(algorithms == BruteAggreg(convScores$ranks, subspace.size,
                                     distance=distance)$top.list[1])]]
        else if(nvm > 1 & subspace.size > 8)
            if(weighted)
                fittedModels[[k]] <- Res[[which(algorithms == RankAggreg(convScores$ranks,
                                     subspace.size, convScores$weights, distance=distance, verbose=FALSE)$top.list[1])]]
            else
                fittedModels[[k]] <- Res[[which(algorithms == RankAggreg(convScores$ranks, subspace.size,
                                     distance=distance, verbose=FALSE)$top.list[1])]]        
        else
            fittedModels[[k]] <- Res[[which.max(scores[, 1])]]
	 # Variable importance evaluation
		if(predictor.importance){
                   model.final <- fittedModels[[k]]$glmnet.fit
                   LassoM.coef <- coef(model.final, s=fittedModels[[k]]$lambda.min)
                   Var_subset <- names(LassoM.coef[as.vector(LassoM.coef[,1]!=0), ])
                   Var_subset <- Var_subset[Var_subset!=c("(Intercept)")]
                   if(!is.null(Var_subset)){
                    training_2 <- training[,Var_subset]
                    Res1 <- cv.glmnet(x=as.matrix(training_2), y=as.factor(trainY), type.measure="deviance", nfolds=kfold, family="binomial")
                    model.final1 <- Res1$glmnet.fit
                    LassoM.coef1 <- coef(model.final1, s=Res1$lambda.min)[-1,]
                        if (length(LassoM.coef1)!=0){
                        for(i0 in 1:length(LassoM.coef1)){
                          for(j0 in 1:length(colnames(RecordM))){
                            if (names(LassoM.coef1)[i0]==colnames(RecordM)[j0]){RecordM[k, j0]=LassoM.coef1[i0]}
                            else {RecordM[k, j0]=RecordM[k, j0]}
                                                                 }
                                                         }
                                                    }
                                           }
		          }
        # Trimmed bagging
		if(trimmed){
                glmFit <- fittedModels[[k]]
                model.glmFit <- glmFit$glmnet.fit
                model.var <- rownames(as.matrix(coef(model.glmFit, s=glmFit$lambda.min)))
                model.testing <- as.matrix(testing[, model.var[-1]]) 
                model.predicted <- predict(glmFit, newx=model.testing, s=c("lambda.min"), type=c("class"))
                model.accuracy[k] <- accuracy(truth, model.predicted)
		           }
	 # Running message
        if(verbose)
            cat("Iter ", k, "\n")
    } # Loop End 1:M
    # Variable importance socres
    RecordM <- abs(RecordM)
    varImportance <- as.matrix(apply(RecordM, 2, mean), ncol(RecordM), 1)
    # Trimmed bagging
    if(trimmed){
        trimmedModels <- fittedModels
          for(r in 1:M){
            trimmedModels[[rank(model.accuracy)[r]]] <- fittedModels[[r]]
		           }
               }
                           } # logistic regression model end

    result <- list(family=family, M=M, predictor.subset=predictor.subset, subspace.size=subspace.size, validation.metric=validation, boot.scale=boot.scale, 
                   distance=distance, models.fitted=fittedModels, models.trimmed=trimmedModels, y.true=y, conv.scores=convScores, importance=varImportance)
    class(result) <- "bagging"
    result
}