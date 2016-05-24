#' Fit linear discriminant
#'
#' Wrapper for the \pkg{MASS} package implementation.
#'
#' @param x Dataset, numerical matrix with observations as rows.
#' @param y Class labels, factor.
#' @param ... Sent to \code{\link{lda}}.
#' @return Fitted linear discriminant.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{predict_lda}},
#'   \code{\link{modeling_procedure}}
#' @export
fit_lda <- function(x, y, ...) {
    nice_require("MASS")
    MASS::lda(x, y, ...)
}

#' Prediction using already trained prediction model
#'
#' Wrapper for the \pkg{MASS} package implementation.
#'
#' @param object Fitted classifier as produced by \code{\link{evaluate}}.
#' @param x Dataset of observations to be classified.
#' @param ... Sent to \code{\link{predict.lda}}.
#' @return A list with elements:
#' \itemize{
#'     \item{\code{prediction}: Factor of predicted class memberships.}
#'     \item{\code{probability}: Data frame of predicted class probabilities.}
#' }
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_lda}},
#'   \code{\link{modeling_procedure}}
#' @export
predict_lda <- function(object, x, ...){
    nice_require("MASS")
    prediction <- predict(object, newdata=x, ...)
    return(list(prediction = prediction$class,
                probability = as.data.frame(prediction$posterior)))
}

