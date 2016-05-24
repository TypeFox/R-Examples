#' Fit a decision tree
#' 
#' @param x Data set (features).
#' @param y Response.
#' @param ... Sent to \code{\link{rpart}}.
#' @return A fitted decision tree.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
fit_rpart <- function(x, y, ...){
    nice_require("rpart", "is needed to fit decision trees")
    model <- if(inherits(y, "formula")){
        rpart::rpart(formula = y, data = x, ...)
    } else {
        rpart::rpart(formula = y ~ ., data = x, ...)
    }
    model$y <- y
    model
}

#' Predict using a fitted decision tree
#' 
#' @param object Fitted decision tree.
#' @param x New data whose response is to be predicted.
#' @return Predictions. The exact form depends on the type of application
#'   (classification or regression)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
predict_rpart <- function(object, x){
    if(is.factor(object$y)){
        # Classification
        list(prediction = predict(object, x, type="class"),
             probability = as.data.frame(predict(object, x, type="prob")))
    } else {
        # Regression
        list(prediction = predict(object, x, type="vector"))
    }
}

