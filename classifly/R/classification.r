#' Generate classification data.
#'
#' Given a model, this function generates points within
#' the range of the data, classifies them, and attempts to locate boundaries
#' by looking at advantage.
#'
#' If posterior probabilities of classification are available, then the
#' \code{\link{advantage}} will be calculated directly.  If not,
#' \code{\link{knn}} is used calculate the advantage based on the number of
#' neighbouring points that share the same classification.  Because knn is
#' $O(n^2)$ this method is rather slow for large (>20,000 say) data sets.
#'
#' By default, the boundary points are identified
#' as those below the 5th-percentile for advantage.
#'
#' @param model classification model
#' @param data data set used in model
#' @param n number of points to generate
#' @param method method to use, currently either grid (an evenly spaced grid),
#'   random (uniform random distribution across cube), or nonaligned (grid +
#'   some random peturbationb)
#' @param advantage if \code{TRUE}, compute advantage, otherwise don't
#' @return data.frame of classified data
#' @keywords datagen
#' @export
generate_classification_data <- function(model, data, n, method, advantage) {
	v <- variables(model)

	df <- generate_data(data[, v$predictors, drop=FALSE], n=n, method=method)
	post <- posterior(model, df)

	df[[".ADVANTAGE"]] <- NA
	if (is.null(post)) {
		df[[v$response]] <- classify(model, df)
		if (advantage) {
			v <- variables(model)
			pred <- rescaler(df[, v$predictors], type="range")
			a <- class::knn(pred, pred, df[,v$response], prob=T, k=5)

			df[[".ADVANTAGE"]] <- attr(a, "prob")

		}
	} else {
		df[[v$response]] <- factor(max.col(post), levels=1:ncol(post), labels=colnames(post))
		if (advantage) df[[".ADVANTAGE"]] <- advantage(post)
		df <- cbind(df, post)
	}
	df[[".TYPE"]] <- factor("simulated")

	df
}

#' Extract classifications from a variety of methods.
#'
#' If the classification method can produce a matrix of posterior
#' probabilities (see \code{\link{posterior}}), then that will be used to
#' calculate the \code{\link{advantage}}.  Otherwise, the classify method
#' will be used and the advantage calculated using a k-nearest neighbours
#' approach.
#'
#' @param model model object
#' @param data data set used in model
#' @param ... other argument passed on to methods
#' @keywords internal
#' @export
classify <- function(model, data, ...) UseMethod("classify", model)
#' @export
classify.rpart <- function(model, data, ...) {
  stats::predict(model, data, type="class")
}


#' Extract posterior group probabilities
#'
#' Every classification method seems to provide a slighly different
#' way of retrieving the posterior probability of group membership.  This
#' function provides a common interface to all of them
#'
#' @param model model object
#' @param data data set used in model
#' @export
posterior <- function(model, data) UseMethod("posterior", model)
#' @export
posterior.default <- function(model, data) NULL
#' @export
posterior.lda <- function(model, data) {
  stats::predict(model, data)$posterior
}
#' @export
posterior.qda <- posterior.lda
#' @export
posterior.randomForest <- function(model, data) {
  stats::predict(model, data, type="prob")
}
#' @export
posterior.svm <- function(model, data) {
  attr(stats::predict(model, data, probability = TRUE), "probabilities")
}
#' @export
posterior.nnet <- function(model, data) {
	probs <- stats::predict(model, data)
	cbind(probs, 1 - rowSums(probs))
}
#' @export
posterior.glm <- function(model, data) {
	probs <- stats::predict(model, data, type="response")
	probs <- cbind(probs, 1 - probs)
	colnames(probs) <- levels(model$model[[variables(model)$response]])
	probs
}

#' Calculate the advantage the most likely class has over the next most
#' likely.
#'
#' This is used to identify the boundaries between classification regions.
#' Points with low (close to 0) advantage are likely to be near boundaries.
#'
#' @param post matrix of posterior probabilities
#' @keywords classif
#' @export
advantage <- function(post) {
	apply(post, 1, function(x) -diff(x[order(x, decreasing=TRUE)[1:2]]))
}

#' Extract predictor and response variables for a model object.
#'
#' Due to the way that most model objects are stored, you
#' also need to supply the data set you used with the original
#' data set.  It currently doesn't support models fitted without
#' using a data argument.
#'
#' @export
#' @param model model object
#' @return list containing response and predictor variables
#' @keywords attribute
variables <- function(model) UseMethod("variables", model)
#' @export
variables.default <- function(model) {
	list(
		response = all.vars(model$terms[[2]]),
		predictors = all.vars(model$terms[[3]])
	)
}
