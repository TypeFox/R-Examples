#' Calculate the Neyman-Pearson Classifier from a sample of class 0 and class 1.
#'
#' \code{npc} calculate the Neyman-Pearson Classifier for a given type I error
#' constraint.
#' @export
#' @importFrom e1071 svm
#' @importFrom e1071 naiveBayes
#' @importFrom glmnet cv.glmnet
#' @importFrom MASS lda
#' @importFrom randomForest randomForest
#' @importFrom ada ada
#' @importFrom parallel mcmapply
#' @importFrom graphics plot
#' @importFrom stats approx
#' @importFrom stats glm
#' @importFrom stats predict
#' @param x n * p observation matrix. n observations, p covariates.
#' @param y n 0/1 observatons.
#' @param method classification method.
#' \itemize{
#' \item logistic: \link{glm} function with family = 'binomial'
#' \item penlog: \code{\link[glmnet]{glmnet}} in \code{glmnet} package
#' \item svm: \code{\link[e1071]{svm}} in \code{e1071} package
#' \item randomforest: \code{\link[randomForest]{randomForest}} in \code{randomForest} package
#' \item lda: \code{\link[MASS]{lda}} in \code{MASS} package
#' \item nb: \code{\link[e1071]{naiveBayes}} in \code{e1071} package
#' \item ada: \code{\link[ada]{ada}} in \code{ada} package
#' \item custom: a custom classifier. score vector needed.
#' }
#' @param kernel kernel used in the svm method. Default = 'radial'.
#' @param score score vector corresponding to y. Required when method  = 'custom'.
#' @param pred.score predicted score vector for the test sample. Optional when method  = 'custom'.
#' @param alpha the desirable control on type I error. Default = 0.05.
#' @param delta the violation rate of the type I error. Default = 0.05.
#' @param split the number of splits for the class 0 sample. Default = 1. For ensemble
#' version, choose split > 1.  When method = 'custom',  split = 0 always.
#' @param loc.prob the precalculated threshold locations in probability.
#'   Default = NULL.
#' @param n.cores number of cores used for parallel computing. Default = 1.
#' @param randSeed the random seed used in the algorithm.
#' @return An object with S3 class npc.
#' \item{fit}{the fit from the specified classifier.}
#' \item{score}{the score vector for each observation.}
#'  \item{cutoff}{thecutoff determined via bootstrap to achieve the specified type I error
#'   control.}
#'  \item{sign}{whether class 1 has a larger average score than the
#'   class 0.}
#'  \item{method}{the classification method.}
#'  \item{loc.prob}{the
#'   percentile used to determine the cutoff for the specified type I error
#'   control.}
#' @seealso \code{\link{nproc}} and \code{\link{predict.npc}}
#' @examples
#' set.seed(0)
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' xtest = matrix(rnorm(n*2),n,2)
#' ctest = 1+3*xtest[,1]
#' ytest = rbinom(n,1,1/(1+exp(-ctest)))
#'
#' ##Use svm classifier and the default type I error control with alpha=0.05
#' fit = npc(x, y, method = 'svm')
#' pred = predict(fit,xtest)
#' fit.score = predict(fit,x)
#' accuracy = mean(pred$pred.label==ytest)
#' cat('Overall Accuracy: ',  accuracy,'\n')
#' ind0 = which(ytest==0)
#' typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' cat('Type I error: ', typeI, '\n')
#'
#' ##Ensembled svm classifier with split = 11,  alpha=0.05
#' #fit = npc(x, y, method = 'svm', split = 11)
#' #pred = predict(fit,xtest)
#' #fit.score = predict(fit,x)
#' #accuracy = mean(pred$pred.label==ytest)
#' #cat('Overall Accuracy: ',  accuracy,'\n')
#' #ind0 = which(ytest==0)
#' #typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' #cat('Type I error: ', typeI, '\n')
#'
#' ##Now, change the method to logistic regression and change alpha to 0.1
#' #fit = npc(x, y, method = 'logistic', alpha = 0.1)
#' #pred = predict(fit,xtest)
#' #accuracy = mean(pred$pred.label==ytest)
#' #cat('Overall Accuracy: ',  accuracy,'\n')
#' #ind0 = which(ytest==0)
#' #typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' #cat('Type I error: ', typeI, '\n')
#'
#' ##Now, change the method to adaboost
#' #fit = npc(x, y, method = 'ada', alpha = 0.1)
#' #pred = predict(fit,xtest)
#' #accuracy = mean(pred$pred.label==ytest)
#' #cat('Overall Accuracy: ',  accuracy,'\n')
#' #ind0 = which(ytest==0)
#' #typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' #cat('Type I error: ', typeI, '\n')
#'
#' ##A 'custom' npc classifier with y and score.
#' #fit2 = npc(y = y, score = fit.score$pred.score,
#' #pred.score = pred$pred.score, method = 'custom')

npc <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest",
                "lda", "nb", "ada", "custom"), kernel = "radial", score = NULL, pred.score = NULL,
                alpha = 0.05, delta = 0.05, split = 1, loc.prob = NULL,
                n.cores = 1, randSeed = 0) {
  method = match.arg(method)
  set.seed(randSeed)
  if (method == "custom" & is.null(score)) {
    stop("score is needed when specifying method = \"custom\"")
  }
  if (!is.null(score)) { ##custom method, user specifed score vector
    test.score = score
    test.y = y
    pred.score = pred.score
    obj = npc.core(test.y, test.score, pred.score = pred.score, alpha = alpha, delta = delta,
                   loc.prob = loc.prob, n.cores = n.cores)
    object = list(pred.y = obj$pred.y, y = test.y, score = test.score,
                  cutoff = obj$cutoff, sign = obj$sign, method = method, loc.prob = loc.prob)
    class(object) = "npc"
    return(object)
   } else {
    n = nrow(x)
    p = ncol(x)
    ind0 = which(y == 0)  ##indices for class 0
    ind1 = which(y == 1)  ##indices for class 1
    n0 = length(ind0)
    n1 = length(ind1)
    object.all = NULL
      for(i in 1:max(1,split)){
        if(split==0){
          ind01 = ind0
          ind02 = ind0
        } else {
          ind01 = sample(ind0, round(n0)/2)  ##use for train classifer
          ind02 = setdiff(ind0, ind01)  ##use for determine cutoff
        }
        train.ind = c(ind01, ind1)
        test.ind = c(ind02, ind1)
        train.x = as.matrix(x[train.ind, ])
        train.y = y[train.ind]
        test.x = as.matrix(x[test.ind, ])
        test.y = y[test.ind]
        colnames(train.x) = paste("x", 1:p, sep = "")
        colnames(test.x) = paste("x", 1:p, sep = "")
        if (method == "logistic") {
          train.data = data.frame(train.x, y = train.y)
          fit = glm(y ~ ., data = train.data, family = "binomial")
          test.score = predict(fit, data.frame(test.x), type = "response")
        } else if (method == "penlog") {
          fit = cv.glmnet(train.x, train.y, family = "binomial")
          test.score = predict(fit$glmnet.fit, newx = test.x, type = "response",
                               s = fit$lambda.min)
          test.score = as.vector(test.score)

        } else if (method == "svm") {
          train.y = as.factor(train.y)
          fit = svm(train.x, train.y, kernel = kernel)
          test.score = attr(predict(fit, test.x, decision.values = TRUE), "decision.values")[,
                                                                                             1]
        } else if (method == "randomforest") {
          train.y = as.factor(train.y)
          fit = randomForest(train.x, train.y)
          test.score = predict(fit, test.x, type = "prob")[, 2]
        } else if (method == "lda") {
          fit = lda(train.x, train.y)
          test.score = predict(fit, data.frame(test.x))$posterior[, 2]
        } else if (method == "nb") {
          fit = naiveBayes(train.x, train.y)
          test.score = predict(fit, data.frame(test.x), type = "raw")[, 2]
        } else if (method == "ada") {
          train.data = data.frame(train.x, y = train.y)
          fit = ada(y~., data=train.data)
          test.score = predict(fit, data.frame(test.x), type = "probs")[, 2]
        }




      obj = npc.core(test.y, test.score, pred.score = pred.score, alpha = alpha, delta = delta,
                     loc.prob = loc.prob, n.cores = n.cores)
      if (is.null(loc.prob)){
        loc.prob = obj$loc.prob
      }
      object = list(fit = fit, pred.y = obj$pred.y, y = test.y, score = test.score,
                    cutoff = obj$cutoff, sign = obj$sign, method = method, loc.prob = loc.prob)
      object.all[[i]] = object
      }
      if(split>1){
      class(object.all) = "npc"
      object.all$split = split
      object.all$method = method
      return(object.all)
      } else{
        class(object) = "npc"
        object$split = split
        return(object)
      }
    }


}
