#' Fit a support vector machine
#' 
#' @param x Dataset, observations as rows.
#' @param y Response vector.
#' @param probability If \code{FALSE} support for class probability estimation
#'   is not included in the fitted model. This may save some computation time.
#' @param ranges Parameter ranges to tune over. Sent to
#'   \code{\link[e1071]{tune.svm}}.
#' @param ... Sent to \code{\link[e1071]{svm}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
fit_svm <- function(x, y, probability=TRUE, ranges, ...){
    nice_require("e1071")
    if(missing(ranges)){
        e1071::svm(x, y, probability=probability, ...)
    } else {
        e1071::tune.svm(x, y, ranges=ranges, probability=probability, ...)
    }
}

#' Predict using support vector machine
#' 
#' @param object Fitted SVM, as produced by \code{\link{fit_svm}}.
#' @param x Data set to predict response for.
#' @param probability Whether to calculate class probabilities.
#' @param statistic Whether to return the raw classification statistics.
#' @param ... Sent to \code{\link[e1071:predict.svm]{predict}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
predict_svm <- function(object, x, probability=object$compprob, statistic=TRUE, ...){
    prediction <- list(prediction, predict(object, x, probability=probability,
                                           decision.values = statistic, ...))
    if(probability){
        prediction$probability <- attr(prediction$prediction, "probabilities")
        attr(prediction$prediction, "probability") <- NULL
    }
    if(statistic){
        prediction$statistic <- attr(prediction$prediction, "decision.value")
        attr(prediction$prediction, "decision.value") <- NULL
    }
    prediction
}

#' Fit a naive Bayes classifier
#' 
#' @param x Dataset, observations as rows.
#' @param y Response vector.
#' @param ... Send to \code{\link[e1071]{naiveBayes}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
fit_naive_bayes <- function(x, y, ...){
    nice_require("e1071")
    object <- e1071::naiveBayes(x, y, ...)
}

#' Predict using naive Bayes model
#' 
#' @param object Fitted naive Bayes model.
#' @param x Data set to predict response for
#' @seealso \code{\link[e1071:predict.svm]{predict}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
predict_naive_bayes <- function(object, x){
    list(prediction = predict(object=object, newdata=x),
         probability = predict(object=object, newdata=x, type="raw"))
}

