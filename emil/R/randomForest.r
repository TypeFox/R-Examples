#' Fit random forest.
#'
#' Directly calling the \pkg{randomForest} package implementation. See
#' \code{\link[randomForest]{randomForest}} for parameter specification.
#'
#' @param x Dataset, numerical matrix with observations as rows.
#' @param y Class labels, factor.
#' @param ... Sent to \code{\link[randomForest]{randomForest}}.
#' @return Fitted random forest.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{predict_randomForest}},
#'   \code{\link{importance_randomForest}}, \code{\link{modeling_procedure}}
#' @export
fit_randomForest <- function(x, y, ...){
    nice_require("randomForest")
    tryCatch(na.fail(x), error = function(...){
        stop("Random forest does not accept any missing values.")
    })
    randomForest::randomForest(x, y, ...)
}


#' Prediction using random forest.
#'
#' @param object Fitted model.
#' @param x Dataset of observations to be classified.
#' @param ... Ignored
#' @return When used for classification, a list with elements:
#' \itemize{
#'     \item{\code{prediction}: Factor of predicted class memberships.}
#'     \item{\code{probability}: Data frame of predicted class probabilities.}
#' }
#'
#' When used for regression, a list with the element:
#' \itemize{
#'     \item{\code{prediction}: Vector of predicted response.}
#' }
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_randomForest}},
#'   \code{\link{importance_randomForest}}, \code{\link{modeling_procedure}}
#' @export
predict_randomForest <- function(object, x, ...){
    nice_require("randomForest")
    if(is.factor(object$y)){
        list(prediction = predict(object, newdata=x, type="response"),
             probability = as.data.frame(predict(object, newdata=x, type="prob")))
    } else {
        list(prediction = predict(object, newdata=x, type="response"))
    }
}


#' Feature importance of random forest.
#' 
#' @param object Fitted randomForest classifier
#' @param type Importance can be assessed in two ways:
#'   \describe{
#'     \item{1.}{Permuted out-of-bag prediction error (default). This can only be
#'            used if the classifier was fitted with argument
#'            \code{prediction=TRUE} which is default.}
#'     \item{2.}{Total decrease in node impurity.}
#'   }
#' @param ... Ignored.
#' @return An prediction vector with elements corresponding to variables.
#' @seealso \code{\link{emil}}, \code{\link{fit_randomForest}},
#'   \code{\link{predict_randomForest}}, \code{\link{modeling_procedure}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
importance_randomForest <- function(object, type=1, ...){
    if(is_blank(object$importance))
        stop("To calculate feature importance of random forsests you must set the fitting parameter `importance=TRUE`.")
    object$importance
}

