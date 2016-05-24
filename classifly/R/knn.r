#' A wrapper function for \code{\link[class]{knn}} to allow use
#' with classifly.
#'
#' @param formula classification formula
#' @param data training data set
#' @param k number of neighbours to use
#' @keywords classif
#' @export
knnf <- function(formula, data, k=2) {
	structure(list(terms=terms(formula, data=data), data=data, k=k),
	   class="knnf")
}

#' @export
classify.knnf <- function(model, data, ...) {
	v <- variables(model)
	class::knn(model$data[,v$predictors], data[,v$predictors],
    model$data[, v$response], k = model$k)
}
