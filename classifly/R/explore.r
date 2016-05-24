#' Classifly provides a convenient method to fit a classification function
#' and then explore the results in the original high dimensional space.
#'
#' This is a convenient function to fit a classification function and
#' then explore the results using GGobi.  You can also do this in two
#' separate steps using the classification function and then
#' \code{\link{explore}}.
#'
#' By default in GGobi, points that are not on the boundary (ie. that have an
#' advantage greater than the 5% percentile) are hidden.  To show them, switch
#' to brush mode and choose include shadowed points from the brush menu on
#' the plot window.  You can then brush them yourself to explore how the
#' certainty of classification varies throughout the space
#'
#' Special notes:
#'
#' \itemize{
#' 	\item You should make sure the response variable is a factor
#' 	\item For SVM, make sure to include \code{probability = TRUE} in the
#'       arguments to \code{classifly}
#'
#' }
#'
#' @param data Data set use for classification
#' @param model Classification formula, usually of the form
#'   \code{response ~ predictors}
#' @param classifier Function to use for the classification, eg.
#'    \code{\link[MASS]{lda}}
#' @param ...  Other arguments passed to classification function.  For
#'    example. if you use \code{\link[e1071]{svm}} you need to use
#'    \code{probabiltiy = TRUE} so that posterior probabilities can be
#'    retrieved.
#' @param n Number of points to simulate.  To maintain the illusion of a
#'   filled solid this needs to increase with dimension.  10,000 points seems
#'   adequate for up to four of five dimensions, but if you have more
#'   predictors than that, you will need to increase this number.
#' @param method method to simulate points: grid, random or nonaligned
#'    (default).  See \code{\link{simvar}} for more details on the methods
#'    used.
#' @param type type of scaling to apply to data.  Defaults to commmon range.
#'   See \code{\link[reshape]{rescaler}} for more details.
#' @aliases classifly package-classifly
#' @seealso \code{\link{explore}}, \url{http://had.co.nz/classifly}
#' @keywords dynamic
#' @export
#' @examples
#' data(kyphosis, package = "rpart")
#' library(MASS)
#' classifly(kyphosis, Kyphosis ~ . , lda)
#' classifly(kyphosis, Kyphosis ~ . , qda)
#' classifly(kyphosis, Kyphosis ~ . , glm, family="binomial")
#' classifly(kyphosis, Kyphosis ~ . , knnf, k=3)
#'
#' library(rpart)
#' classifly(kyphosis, Kyphosis ~ . , rpart)
#'
#' \donttest{
#' if (require("e1071")) {
#' classifly(kyphosis, Kyphosis ~ . , svm, probability=TRUE)
#' classifly(kyphosis, Kyphosis ~ . , svm, probability=TRUE, kernel="linear")
#' classifly(kyphosis, Kyphosis ~ . , best.svm, probability=TRUE,
#'    kernel="linear")
#'
#' # Also can use explore directorly
#' bsvm <- best.svm(Species~., data = iris, gamma = 2^(-1:1),
#'   cost = 2^(2:+ 4), probability=TRUE)
#' explore(bsvm, iris)
#' }
#' }
classifly <- function(data, model, classifier, ..., n=10000, method="nonaligned", type="range") {
  data <- rescaler(data, type=type)
	classifly <- classifier(model, data=data, ...)
	explore(classifly, data, n=n, method=method, advantage=TRUE)
}

#' Default method for exploring objects
#'
#' The default method currently works for classification
#' functions.
#'
#' It generates a data set filling the design space, finds
#' class boundaries (if desired) and then displays in a new
#' ggobi instance.
#'
#' @param model classification object
#' @param data data set used with classifier
#' @param n number of points to generate when searching for boundaries
#' @param method method to generate points, see \code{\link{generate_data}}
#' @param advantage only display boundaries
#' @param ... other arguments not currently used
#' @seealso \code{\link{generate_classification_data}},
#'    \url{http://had.co.nz/classifly}
#' @return A \code{\link{invisible}} data frame of class \code{classifly}
#'   that contains all the simulated and true data.  This can be saved and
#'   then printed later to open with rggobi.
#' @export
#' @examples
#' if (require("e1071")) {
#' bsvm <- best.svm(Species~., data = iris, gamma = 2^(-1:1),
#'   cost = 2^(2:+ 4), probability=TRUE)
#' explore(bsvm, iris)
#' }
explore <- function(model, data, n=10000, method="nonaligned", advantage=TRUE, ...) {
	v <- variables(model)
	grid <- generate_classification_data(model, data, n=n, method=method, advantage=TRUE)
	actual <- data[,c(v$predictor, v$response)]
	actual[[".TYPE"]] <- factor("actual")

	data <- plyr::rbind.fill(grid, actual)

	boundary <- stats::quantile(data[[".ADVANTAGE"]], 0.1, na.rm=TRUE)

	data$.BOUNDARY <- !is.na(data[[".ADVANTAGE"]]) & data[[".ADVANTAGE"]] > boundary

	class(data) <- c("classifly", class(data))
	attr(data, "variables") <- v
	attr(data, "boundary") <- boundary
	data
}

#' @export
print.classifly <- function(x, ...) {
  if (!interactive()) {
    message("Skipping interactive exploration")
    return(invisible())
  }

	if (!require("rggobi", quietly=TRUE))
    stop("rggobi required to visualise classifications in GGobi.")

	v <- attr(x, "variables")
	g <- rggobi::ggobi(x)

	d <- g[1]
  rggobi::glyph_colour(d) <- as.numeric(x[[v$response]]) + 1
  rggobi::glyph_type(d) <- ifelse(x[[".TYPE"]] == "simulated", 1, 6)
  rggobi::shadowed(d) <- x$.BOUNDARY
  rggobi::excluded(d) <- x$.BOUNDARY
	invisible(d)
}
