#' Build transformation models
#'
#' Base classifiers are used to build models to solve the the transformation
#' problems. To create a new base classifier, two steps are necessary:
#' \enumerate{
#'   \item Create a train method
#'   \item Create a prediction method
#' }
#' This section is about how to create the first step: a train method.
#' To create a new predict model see \code{\link{mlpredict}} documentation.
#'
#' @section How to create a new train base method:
#' First, is necessary to define a name of your classifier, because this name
#' determines the method name. The base method name must start with
#' \code{mltrain.base} followed by the designed name, e.g. a \code{'FOO'}
#' classify must be defined as \code{mltrain.baseFOO} (we suggest always use
#' upper case names).
#'
#' Next, your method must receive at least two parameters (\code{object, ...}).
#' Use \code{object$data[, object$labelindex]} or
#' \code{object$data[, object$labelname]} to access the labels values and use
#' \code{object$data[, -object$labelindex]} to access the predictive attributes.
#' If you need to know which are the multi-label dataset and method, use
#' \code{object$mldataset} and \code{object$mlmethod}, respectively.
#'
#' Finally, your method should return a model that will be used by the mlpredict
#' method. Remember, that your method may be used to buid binary and multi-class
#' models.
#'
#' @param object A \code{mltransformation} object. This is used as a list and
#' contains at least five values:
#'  \describe{
#'    \item{object$data}{A data.frame with the train data, where the columns are
#'    the attributes and the rows are the examples.}
#'    \item{object$labelname}{The name of the class column.}
#'    \item{object$labelindex}{The column index of the class.}
#'    \item{object$mldataset}{The name of multi-label dataset.}
#'    \item{object$mlmethod}{The name of the multi-label method.}
#'  }
#'  Others values may be specified by the multi-label method.
#' @param ... Others arguments passed to the base method.
#' @return A model object. The class of this model can be of any type, however,
#'  this object will be passed to the respective mlpredict method.
#' @export
#'
#' @examples
#' # Create a empty model of type FOO
#' mltrain.baseFOO <- function (object, ...) {
#'    mymodel <- list(
#'      classes = as.character(unique(object$data[, object$labelindex]))
#'    )
#'    class(mymodel) <- 'fooModel'
#'    mymodel
#' }
#'
#' # Using this base method with Binary Relevance
#' brmodel <- br(toyml, 'FOO')
#'
#' \dontrun{
#'
#' # Create a SVM method using the e1071 package
#' library(e1071)
#' mltrain.baseSVM <- function (object, ...) {
#'    traindata <- object$data[, -object$labelindex]
#'    labeldata <- object$data[, object$labelindex]
#'    model <- svm(traindata, labeldata, probability = TRUE, ...)
#'    model
#' }
#' }
mltrain <- function(object, ...) {
  UseMethod("mltrain")
}

#' Prediction transformation problems
#'
#' Base classifiers are used to build models to solve the the transformation
#' problems. To create a new base classifier, two steps are necessary:
#' \enumerate{
#'   \item Create a train method
#'   \item Create a prediction method
#' }
#' This section is about how to create the second step: a prediction method.
#' To create a new train method see \code{\link{mltrain}} documentation.
#'
#' @section How to create a new prediction base method:
#' Fist is necessary to know the class of model generate by the respective train
#' method, because this name determines the method name. It must start with
#' \code{'mlpredict.'}, followed by the model class name, e.g. a model with
#' class 'fooModel' must be called as \code{mlpredict.fooModel}.
#'
#' After defined the name, you need to implement your prediction base method.
#' The model built on mltrain is available on \code{model} parameter and the
#' \code{newdata} is the data to be predict.
#'
#' The return of this method must be a data.frame with two columns called
#' \code{"prediction"} and \code{"probability"}. The first column contains the
#' predicted classe and the second the probability/score/confidence of this
#' prediction. The rows represents the examples.
#'
#' @param model An object model returned by some mltrain method, its class
#'  determine the name of this method.
#' @param newdata A data.frame with the new data to be predicted.
#' @param ... Others arguments passed to the predict method.
#' @return A matrix with the probabilities of each class value/example,
#'  where the rows are the examples and the columns the class values.
#' @export
#'
#' @examples
#'
#' # Create a method that predict always the first class
#' # The model must be of the class 'fooModel'
#' mlpredict.fooModel <- function (model, newdata, ...) {
#'    # Predict the first class with a random confidence
#'    data.frame(
#'      prediction = rep(model$classes[1], nrow(newdata)),
#'      probability = sapply(runif(nrow(newdata)), function (score) {
#'        max(score, 1 - score)
#'      }),
#'      row.names = rownames(newdata)
#'    )
#' }
#'
#' \dontrun{
#' # Create a SVM predict method using the e1071 package (the class of SVM model
#' # from e1071 package is 'svm')
#' library(e1071)
#' mlpredict.svm <- function (dataset, newdata, ...) {
#'    result <- predict(model, newdata, probability = TRUE, ...)
#'    attr(result, 'probabilities')
#' }
#' }
mlpredict <- function(model, newdata, ...) {
  UseMethod("mlpredict")
}

# DEFAULT METHOD -------------------------------------------------------------
#' @describeIn mltrain Default S3 method
#' @export
mltrain.default <- function(object, ...) {
  funcname <- paste("mltrain.base", object$methodname, sep = "")
  stop(paste("The function '", funcname, "(object, ...)' is not implemented",
             sep = ""))
}

#' @describeIn mlpredict Default S3 method
#' @export
mlpredict.default <- function(model, newdata, ...) {
  funcname <- paste("mlpredict.", class(model), sep = "")
  stop(paste("The function '", funcname,
             "(dataset, newdata, ...)' is not implemented", sep = ""))
}

# SVM METHOD ------------------------------------------------------------------
#' @describeIn mltrain SVM implementation (require \pkg{e1071} package to use)
#' @export
mltrain.baseSVM <- function(object, ...) {
  if (requireNamespace("e1071", quietly = TRUE)) {
    traindata <- object$data[, -object$labelindex]
    labeldata <- object$data[, object$labelindex]
    model <- e1071::svm(traindata, labeldata, probability = TRUE, ...)
  }
  else {
    stop(paste("There are no installed package 'e1071' to use SVM classifier",
                "as base method"))
  }
  model
}

#' @describeIn mlpredict SVM implementation (require \pkg{e1071} package to use)
#' @export
mlpredict.svm <- function(model, newdata, ...) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop(paste("There are no installed package 'e1071' to use SVM classifier",
               "as base method"))
  }

  result <- stats::predict(model, newdata, probability = TRUE, ...)
  prediction <- as.character(result)
  all.prob <- attr(result, "probabilities")

  data.frame(
    prediction = prediction,
    probability = all.prob[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# J48 METHOD -------------------------------------------------------------
#' @describeIn mltrain J48 implementation (require \pkg{RWeka} package to use)
#' @export
mltrain.baseJ48 <- function(object, ...) {
  #http://r.789695.n4.nabble.com/How-to-save-load-RWeka-models-into-from-a-file-td870876.html
  #https://github.com/s-u/rJava/issues/25
  #inspect the JVM log error file in the second execution
  if (requireNamespace("RWeka", quietly = TRUE) &&
      requireNamespace("rJava", quietly = TRUE)) {
    formula <- stats::as.formula(paste("`", object$labelname, "` ~ .", sep=""))
    model <- RWeka::J48(formula, object$data, ...)
    rJava::.jcache(model$classifier)
  }
  else {
    stop(paste("There are no installed package 'RWeka' and 'rJava' to use J48",
               "classifier as base method"))
  }
  model
}

#' @describeIn mlpredict J48 implementation (require \pkg{RWeka} package to use)
#' @export
mlpredict.J48 <- function(model, newdata, ...) {
  if (!requireNamespace("RWeka", quietly = TRUE)) {
    stop(paste("There are no installed package 'RWeka' to use J48 classifier",
               "as base method"))
  }

  result <- stats::predict(model, newdata, type = "probability", ...)
  prediction <- colnames(result)[apply(result, 1, which.max)]
  data.frame(
    prediction = prediction,
    probability = result[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# C5.0 METHOD ------------------------------------------------------------------
#' @describeIn mltrain C5.0 implementation (require \pkg{C50} package to use)
#' @export
mltrain.baseC5.0 <- function(object, ...) {
  if (requireNamespace("C50", quietly = TRUE)) {
    traindata <- object$data[, -object$labelindex]
    labeldata <- object$data[, object$labelindex]
    model <- C50::C5.0(traindata, labeldata, ...)
  }
  else {
    stop(paste("There are no installed package 'C50' to use C5.0 classifier",
               "as base method"))
  }
  model
}

#' @describeIn mlpredict C5.0 implementation (require \pkg{C50} package to use)
#' @export
mlpredict.C5.0 <- function(model, newdata, ...) {
  if (!requireNamespace("C50", quietly = TRUE)) {
    stop(paste("There are no installed package 'C50' to use C5.0 classifier",
               "as base method"))
  }
  result <- C50::predict.C5.0(model, newdata, type = "prob", ...)
  prediction <- colnames(result)[apply(result, 1, which.max)]
  data.frame(
    prediction = prediction,
    probability = result[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# CART METHOD -----------------------------------------------------------------
#' @describeIn mltrain CART implementation (require \pkg{rpart} package to use)
#' @export
mltrain.baseCART <- function(object, ...) {
  if (requireNamespace("rpart", quietly = TRUE)) {
    formula <- stats::as.formula(paste("`", object$labelname, "` ~ .", sep=""))
    model <- rpart::rpart(formula, object$data, ...)
  }
  else {
    stop(paste("There are no installed package 'rpart' to use Cart classifier",
               "as base method"))
  }
  model
}

#' @describeIn mlpredict CART implementation (require \pkg{rpart} package)
#' @export
mlpredict.rpart <- function(model, newdata, ...) {
  if (!requireNamespace("rpart", quietly = TRUE))  {
    stop(paste("There are no installed package 'rpart' to use Cart classifier",
               "as base method"))
  }
  result <- stats::predict(model, newdata, type = "prob", ...)
  prediction <- colnames(result)[apply(result, 1, which.max)]
  data.frame(
    prediction = prediction,
    probability = result[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# RANDOM FOREST METHOD --------------------------------------------------------
#' @describeIn mltrain Random Forest (RF) implementation (require
#'  \pkg{randomForest} package to use)
#' @export
mltrain.baseRF <- function(object, ...) {
  if (requireNamespace("randomForest", quietly = TRUE)) {
    traindata <- object$data[, -object$labelindex]
    labeldata <- object$data[, object$labelindex]
    model <- randomForest::randomForest(traindata, labeldata, ...)
  }
  else {
    stop(paste("There are no installed package 'randomForest' to use",
               "randomForest classifier as base method"))
  }
  model
}

#' @describeIn mlpredict Random Forest (RF) implementation (require
#'  \pkg{randomForest} package to use)
#' @export
mlpredict.randomForest <- function(model, newdata, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop(paste("There are no installed package 'randomForest' to use",
               "randomForest classifier as base method"))
  }

  result <- stats::predict(model, newdata,
                                                type = "prob", ...)
  prediction <- colnames(result)[apply(result, 1, which.max)]
  data.frame(
    prediction = prediction,
    probability = result[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# NAIVE BAYES METHOD ----------------------------------------------------------
#' @describeIn mltrain Naive Bayes (NB) implementation (require
#'  \pkg{e1071} package to use)
#' @export
mltrain.baseNB <- function(object, ...) {
  if (requireNamespace("e1071", quietly = TRUE)) {
    traindata <- object$data[, -object$labelindex]
    labeldata <- object$data[, object$labelindex]
    model <- e1071::naiveBayes(traindata, labeldata, type = "raw", ...)
  }
  else {
    stop(paste("There are no installed package 'e1071' to use naiveBayes",
               "classifier as base method"))
  }
  model
}

#' @describeIn mlpredict Naive Bayes (NB) implementation (require
#'  \pkg{e1071} package to use)
#' @export
mlpredict.naiveBayes <- function(model, newdata, ...) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop(paste("There are no installed package 'e1071' to use naiveBayes",
               "classifier as base method"))
  }
  result <- stats::predict(model, newdata, type = "raw", ...)
  rownames(result) <- rownames(newdata)
  prediction <- colnames(result)[apply(result, 1, which.max)]
  data.frame(
    prediction = prediction,
    probability = result[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# KNN METHOD ------------------------------------------------------------------
#' @describeIn mltrain kNN implementation (require \pkg{kknn} package to use)
#' @export
mltrain.baseKNN <- function(object, ...) {
  if (!requireNamespace("kknn", quietly = TRUE)) {
    stop(paste("There are no installed package 'kknn' to use kNN classifier as",
               "base method"))
  }

  object$extrakNN <- list(...)
  object
}

#' @describeIn mlpredict kNN implementation (require \pkg{kknn} package to use)
#' @export
mlpredict.baseKNN <- function(model, newdata, ...) {
  if (!requireNamespace("kknn", quietly = TRUE)) {
    stop(paste("There are no installed package 'kknn' to use kNN classifier as",
               "base method"))
  }

  formula <- stats::as.formula(paste("`", model$labelname, "` ~ .", sep = ""))
  args <- list(...)
  if (is.null(model$extrakNN[["k"]]) || !is.null(args[["k"]])) {
    result <- kknn::kknn(formula, model$data, newdata, ...)
  }
  else {
    result <- kknn::kknn(formula, model$data, newdata,
                         k = model$extrakNN[["k"]], ...)
  }

  prediction <- as.character(result$fitted.values)
  all.prob <- as.matrix(result$prob)
  rownames(all.prob) <- rownames(newdata)
  data.frame(
    prediction = prediction,
    probability = all.prob[cbind(rownames(newdata), prediction)],
    row.names = rownames(newdata)
  )
}

# Majority METHOD ------------------------------------------------------------
#' @describeIn mltrain Majority model
#' @export
mltrain.baseMAJORITY <- function(object, ...) {
  values <- table(object$data[, object$labelindex])
  model <- list(
    classes = names(values),
    predict = names(which.max(values))
  )
  class(model) <- 'majorityModel'
  model
}

#' @describeIn mlpredict Majority prediction
#' @export
mlpredict.majorityModel <- function(model, newdata, ...) {
  data.frame(
    prediction = rep(model$predict, nrow(newdata)),
    probability = rep(1, nrow(newdata)),
    row.names = rownames(newdata)
  )
}

# Random METHOD ------------------------------------------------------------
#' @describeIn mltrain Random model
#' @export
mltrain.baseRANDOM <- function(object, ...) {
  model <- list(
    classes = as.character(unique(object$data[, object$labelindex]))
  )
  class(model) <- 'randomModel'
  model
}

#' @describeIn mlpredict Majority prediction
#' @export
mlpredict.randomModel <- function(model, newdata, ...) {
  data.frame(
    prediction = sample(model$classes, nrow(newdata), replace = TRUE),
    probability = sapply(stats::runif(nrow(newdata)), function (score) {
      max(score, 1 - score)
    }),
    row.names = rownames(newdata)
  )
}

#' Print Majority model
#' @param x The base model
#' @param ... ignored
#' @export
print.majorityModel <- function (x, ...) {
  cat("Majority Base Model\n\n")
  cat("Label: ", attr(x, "label"), "\n")
  cat("Classes: ", paste(x$classes, collapse = ' | '))
  cat("Predict: ", x$predict)
}

#' Print Random model
#' @param x The base model
#' @param ... ignored
#' @export
print.randomModel <- function (x, ...) {
  cat("Random Base Model\n\n")
  cat("Label: ", attr(x, "label"), "\n")
  cat("Classes: ", paste(x$classes, collapse = ' | '))
}
