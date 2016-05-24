##############################################################
# Function to compare some methods for discriminant analysis #
##############################################################

# With this function you can compare some methods for discriminant analysis.
# In each of the ntimes (see below) iterations the given dataset will be
# divided randomly into a training set and a validation set to fit and
# validate the models.

# Usage:
########

# disc.comp(y, x, methods = c("lda", "qda", "logistic", "rf", "boost",
#                   "neural", "sv", "knn"), nfold = 5, ntimes, k.knn = 1)

# Arguments:
############

#       y: binary response variable
#       x: matrix (dimension n times p) of predictor variables
# methods: vector of the methods to be computed; possible values are
#          "lda" (linear discriminant analysis), "qda" (quadratic
#          discriminant analysis), "logistic" (penalized logistic
#          regression), "rf" (random forest), "boost" (Boosting),
#          "neural" (neural networks), "sv" (support vector machine),
#          "knn" (k nearest neighbour classification)
#   nfold: The dataset will be divided randomly into nfold subsets.
#          One of these subsets will be used to validate the fit and
#          the other will be used to fit the models.
#  ntimes: ntimes iterations will be computed
#   k.knn: the parameter "k" for the computation of k nearest neighbour
#          classification

# Value:
########

# list: one entry for every method of "methods"; In each of these
# entries you can find a list of length "ntimes" with the following
# entries: the misclassification rate ("misclass"), the predicted values for
# the response variable of the test data ("predict"), for almost all methods
# the fits ("fit"), for almost all methods the fitted probabilities ("prob") and
# for almost all methods the squared error of prediction (sum(y - pi)) ("sqerr").
# There are three additional entries: "train": a list of length
# ntimes; In entry i you can find the indices of the observations
# used as training data in the iteration i (i = 1,...,ntimes).
# "x.data": "x" and "y" as dataframe; the values of "y" are possibly
# changed to numbers (1, 2, ...). "ntimes": See ntimes above.

# Required packages:
####################

# "MASS", "glmnet", "nnet", "mboost", "party", "class", "e1071"

disc.comp <- function (y, x, methods = c("lda", "qda", "logistic", "rf", "boost", "neural", "sv", "knn"), nfold = 5, ntimes, k.knn = 1) {

  im <- 1
  nameth <- c()
  res <- list()

  y.mat <- matrix(0, nrow = length(y), ncol = length(levels(y)))
  for (i in 1:length(y)) {
    y.mat[i, as.numeric(y[i])] <- 1
  }

  mult <- F
  if (length(levels(y)) > 2) mult <- T

  if (!mult) y.mat <- y.mat[, -1]

  x.data <- data.frame(cbind(x,y))
  names(x.data) <- c(names(x.data[, -ncol(x.data)]), "y")
  x.data$y <- as.factor(as.numeric(x.data$y))
  form <- as.formula(paste("y ~ ", paste(names(x.data)[1:ncol(x)], collapse = "+", sep = ""), sep = ""))

  test.length <- round(length(y) - length(y) / nfold)
  
  # Defining the training datasets:
  train <- list()
  for (i in 1:ntimes) train[[i]] <- sample(1:length(y), test.length)

  if ("lda" %in% methods) {

    require(MASS)

    fit.lda <- list()
    predict.lda <- list()
    misclass.lda <- c()
    prob.lda <- list()
    sqerr.lda <- c()

    for (i in 1:ntimes) {

      fit.lda[[i]] <- lda(form, data = x.data, subset = train[[i]])
      prediction.lda <- predict(fit.lda[[i]], x.data[-train[[i]], ])
      predict.lda[[i]] <- prediction.lda$class
      prob.lda[[i]] <- if (mult) prediction.lda$posterior else prediction.lda$posterior[,2]
      misclass.lda <- c(misclass.lda, sum(predict.lda[[i]] != x.data$y[-train[[i]]]) / test.length)
      sqerr.lda <- if (mult) c(sqerr.lda, sum((y.mat[-train[[i]], ] - prob.lda[[i]])^2)) else c(sqerr.lda, sum((y.mat[-train[[i]]] - prob.lda[[i]])^2))
    }

    res[[im]] <- list("fit" = fit.lda, "predict" = predict.lda, "misclass" = misclass.lda, "prob" = prob.lda, "sqerr" = sqerr.lda)
    nameth <- c(nameth, "lda")
    im <- im + 1
  }

  if ("qda" %in% methods) {

    require(MASS)

    fit.qda <- list()
    predict.qda <- list()
    misclass.qda <- c()
    prob.qda <- list()
    sqerr.qda <- c()

    for (i in 1:ntimes) {

      fit.qda[[i]] <- qda(form, data = x.data, subset = train[[i]])
      prediction.qda <- predict(fit.qda[[i]], x.data[-train[[i]], ])
      predict.qda[[i]] <- prediction.qda$class
      misclass.qda <- c(misclass.qda, sum(predict.qda[[i]] != x.data$y[-train[[i]]]) / test.length)
      prob.qda[[i]] <- if (mult) prediction.qda$posterior else prediction.qda$posterior[,2]
      sqerr.qda <- if (mult) c(sqerr.qda, sum((y.mat[-train[[i]], ] - prob.qda[[i]])^2)) else c(sqerr.qda, sum((y.mat[-train[[i]]] - prob.qda[[i]])^2))
    }

    res[[im]] <- list("fit" = fit.qda, "predict" = predict.qda, "misclass" = misclass.qda, "prob" = prob.qda, "sqerr" = sqerr.qda)
    nameth <- c(nameth, "qda")
    im <- im + 1
  }

  if ("logistic" %in% methods) {

    require(glmnet)

    fit.logistic <- list()
    cv.logistic <- list()
    predict.logistic <- list()
    misclass.logistic <- c()
    prob.logistic <- list()
    sqerr.logistic <- c()

    for (i in 1:ntimes) {
      # A 5-fold cross-validation is done for choosing the optimal lambda:
      cv.logistic[[i]] <- cv.glmnet(x[train[[i]], ], y[train[[i]]], family = if (mult) "multinomial" else "binomial", nfolds = 5)$lambda.min
#      glmnet.fit.seq <- glmnet(x[train[[i]], ], y[train[[i]]], family = if (mult) "multinomial" else "binomial")
#      lam.seq <- glmnet.fit.seq$lambda
#      dev.seq <- glmnet.fit.seq$dev
#      cv.logistic[[i]] <- lam.seq[which(dev.seq == max(dev.seq))]
      
      fit.logistic[[i]] <- glmnet(x[train[[i]], ], y[train[[i]]], family = if (mult) "multinomial" else "binomial", lambda = cv.logistic[[i]])
      predict.logistic[[i]] <- predict(fit.logistic[[i]], x[-train[[i]], ], type = "class", s = cv.logistic[[i]])
      misclass.logistic <- c(misclass.logistic, sum(predict.logistic[[i]] != as.numeric(y[-train[[i]]])) / test.length)
      prob.logistic[[i]] <- if (mult) predict(fit.logistic[[i]], x[-train[[i]], ], type = "response", s = cv.logistic[[i]])[, , 1] else predict(fit.logistic[[i]], x[-train[[i]], ], type = "response", s = cv.logistic[[i]])
      sqerr.logistic <- if (mult) c(sqerr.logistic, sum((y.mat[-train[[i]], ] - prob.logistic[[i]])^2)) else c(sqerr.logistic, sum((y.mat[-train[[i]]] - prob.logistic[[i]])^2))
    }
    
    res[[im]] <- list("fit" = fit.logistic, "predict" = predict.logistic, "misclass" = misclass.logistic, "prob" = prob.logistic, "sqerr" = sqerr.logistic)
    nameth <- c(nameth, "logistic")
    im <- im + 1
  }

  if ("rf" %in% methods) {

    require(party)

    fit.rf <- list()
    predict.rf <- list()
    misclass.rf <- c()

    for (i in 1:ntimes) {
      fit.rf[[i]] <- cforest(form, data = x.data, subset = train[[i]])
      predict.rf[[i]] <- predict(fit.rf[[i]], x.data[-train[[i]], ], OOB = T)
      misclass.rf <- c(misclass.rf, sum(predict.rf[[i]] != x.data$y[-train[[i]]]) / test.length)
    }
    res[[im]] <- list("fit" = fit.rf, "predict" = predict.rf, "misclass" = misclass.rf)
    nameth <- c(nameth, "rf")
    im <- im + 1
  }

  if ("boost" %in% methods & !mult) {

    require(mboost)

    fit.boost <- list()
    predict.boost <- list()
    misclass.boost <- c()
    prob.boost <- list()
    sqerr.boost <- c()

    for (i in 1:ntimes) {
      fit.boost[[i]] <- glmboost(form, data = x.data[train[[i]], ], family = Binomial())
      predict.boost[[i]] <- predict(fit.boost[[i]], x.data[-train[[i]], ], type = "class")
      misclass.boost <- c(misclass.boost, sum(predict.boost[[i]] != x.data$y[-train[[i]]]) / test.length)
      prob.boost[[i]] <- predict(fit.boost[[i]], x.data[-train[[i]], ], type = "response")
      sqerr.boost <- if (mult) c(sqerr.boost, sum((y.mat[-train[[i]], ] - prob.boost[[i]])^2)) else c(sqerr.boost, sum((y.mat[-train[[i]]] - prob.boost[[i]])^2))
    }
    res[[im]] <- list("fit" = fit.boost, "predict" = predict.boost, "misclass" = misclass.boost, "prob" = prob.boost, "sqerr" = sqerr.boost)
    nameth <- c(nameth, "boost")
    im <- im + 1
  }

  if ("neural" %in% methods) {

    require(nnet)

    fit.neural <- list()
    predict.neural <- list()
    misclass.neural <- c()
    prob.neural <- list()
    sqerr.neural <- c()

    for (i in 1:ntimes) {

      fit.neural[[i]] <- multinom(form, x.data, subset = train[[i]])
      predict.neural[[i]] <- predict(fit.neural[[i]], x.data[-train[[i]], ], type = "class")
      misclass.neural <- c(misclass.neural, sum(predict.neural[[i]] != x.data$y[-train[[i]]]) / test.length)
      prob.neural[[i]] <- predict(fit.neural[[i]], x.data[-train[[i]], ], type = "probs")
      sqerr.neural <- if (mult) c(sqerr.neural, sum((y.mat[-train[[i]], ] - prob.neural[[i]])^2)) else c(sqerr.neural, sum((y.mat[-train[[i]]] - prob.neural[[i]])^2))
    }

    res[[im]] <- list("fit" = fit.neural, "predict" = predict.neural, "misclass" = misclass.neural, "prob" = prob.neural, "sqerr" = sqerr.neural)
    nameth <- c(nameth, "neural")
    im <- im + 1
  }

  if ("sv" %in% methods) {

    require(e1071)

    fit.sv <- list()
    predict.sv <- list()
    misclass.sv <- c()
    prob.sv <- list()
    sqerr.sv <- c()

    for (i in 1:ntimes) {

      fit.sv[[i]] <- svm(form, x.data, subset = train[[i]], kernel =
       "radial",probability = T)
      predict.sv[[i]] <- predict(fit.sv[[i]], x.data[-train[[i]], ])
      misclass.sv <- c(misclass.sv, sum(predict.sv[[i]] != x.data$y[-train[[i]]]) / test.length)
      if (mult) {
        prob.sv.data <- as.data.frame(attr(predict(fit.sv[[i]], x.data[-train[[i]], ], probability = T), "probabilities"))
        prob.sv[[i]] <- as.matrix(prob.sv.data[, sort(names(prob.sv.data))])
      } else prob.sv[[i]] <- attr(predict(fit.sv[[i]], x.data[-train[[i]], ], probability = T), "probabilities")[,3 - as.numeric(x.data$y[train[[i]][1]])]
      sqerr.sv <- if (mult) c(sqerr.sv, sum((y.mat[-train[[i]], ] - prob.sv[[i]])^2)) else c(sqerr.sv, sum((y.mat[-train[[i]]] - prob.sv[[i]])^2))
    }

    res[[im]] <- list("fit" = fit.sv, "predict" = predict.sv, "misclass" = misclass.sv, "prob" = prob.sv, "sqerr" = sqerr.sv)
    nameth <- c(nameth, "sv")
    im <- im + 1
  }

  if ("knn" %in% methods) {

    require(class)

    fit.knn <- list()
    predict.knn <- list()
    misclass.knn <- c()

    for (i in 1:ntimes) {
      predict.knn[[i]] <- knn(x[train[[i]], ], x[-train[[i]], ], y[train[[i]]], k.knn)
      misclass.knn <- c(misclass.knn, sum(predict.knn[[i]] != y[-train[[i]]]) / test.length)
    }
    res[[im]] <- list("predict" = predict.knn, "misclass" = misclass.knn)
    nameth <- c(nameth, "knn")
    im <- im + 1
  }

  res[[im]] <- train
  nameth <- c(nameth, "train")
  im <- im + 1
  res[[im]] <- x.data
  nameth <- c(nameth, "x.data")
  im <- im + 1
  res[[im]] <- ntimes
  nameth <- c(nameth, "ntimes")

  names(res) <- nameth
  return(res)
}

####################################
# Function for plotting ROC-curves #
####################################

# This function computes ROC-curves for the results computed by
# "comp.disc()".  Using ordered factors for the response variable,
# the upper one corresponds to the "positives" (see the manual of
# package "ROCR" (Sing, T. et al., 2009, p. 9 - 10), which is
# used in this function, for details).
# Note that the fitted probabilities and validation sets
# of all iterations (ntimes) computed by "comp.disc()" are used
# to construct the ROC-curve(s).

# Arguments:
############

#  disc: an object generated by "disc.comp()"
# roc.k: a vector of names of the methods which
#        should be included in the plot of the ROC-curves
#        (see argument "method" of the function "disc.comp()");
#        Note that you can´t plot the ROC-curves of "rf" and
#        "knn". These methods will be ignored.

# Required package:
###################

# "ROCR"

roc.curve <- function(disc, roc.k) {

  require(ROCR)
  ind <- 1

  if("lda" %in% roc.k) {
    prob <- c()
    y.true <- c()
    for (i in 1:disc$ntimes) {
      prob <- c(prob, disc$lda$prob[[i]])
      y.true <- c(y.true, disc$x.data$y[-disc$train[[i]]])
    }

    pred <- prediction(prob, y.true)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    plot(perf, lty = ind, lwd = 2)
    ind <- ind + 1
  }

  if("qda" %in% roc.k) {
    prob <- c()
    y.true <- c()
    for (i in 1:disc$ntimes) {
      prob <- c(prob, disc$qda$prob[[i]])
      y.true <- c(y.true, disc$x.data$y[-disc$train[[i]]])
    }

    pred <- prediction(prob, y.true)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    plot(perf, lty = ind, add = (ind > 1), lwd = 2)
    ind <- ind + 1
  }

  if("logistic" %in% roc.k) {
    prob <- c()
    y.true <- c()
    for (i in 1:disc$ntimes) {
      prob <- c(prob, disc$logistic$prob[[i]])
      y.true <- c(y.true, disc$x.data$y[-disc$train[[i]]])
    }

    pred <- prediction(prob, y.true)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    plot(perf, lty = ind, add = (ind > 1), lwd = 2)
    ind <- ind + 1
  }

  if("boost" %in% roc.k) {
    prob <- c()
    y.true <- c()
    for (i in 1:disc$ntimes) {
      prob <- c(prob, disc$boost$prob[[i]])
      y.true <- c(y.true, disc$x.data$y[-disc$train[[i]]])
    }

    pred <- prediction(prob, y.true)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    plot(perf, lty = ind, add = (ind > 1), lwd = 2)
    ind <- ind + 1
  }

  if("neural" %in% roc.k) {
    prob <- c()
    y.true <- c()
    for (i in 1:disc$ntimes) {
      prob <- c(prob, disc$neural$prob[[i]])
      y.true <- c(y.true, disc$x.data$y[-disc$train[[i]]])
    }

    pred <- prediction(prob, y.true)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    plot(perf, lty = ind, add = (ind > 1), lwd = 2)
    ind <- ind + 1
  }


  if("sv" %in% roc.k) {
    prob <- c()
    y.true <- c()
    for (i in 1:disc$ntimes) {
      prob <- c(prob, disc$sv$prob[[i]])
      y.true <- c(y.true, disc$x.data$y[-disc$train[[i]]])
    }

    pred <- prediction(prob, y.true)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    plot(perf, lty = ind, add = (ind > 1), lwd = 2)
    ind <- ind + 1
  }
  legend("bottomright", legend = roc.k, lty = 1:ind, lwd = 2)
}
