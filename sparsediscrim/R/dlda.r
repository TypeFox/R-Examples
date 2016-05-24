#' Diagonal Linear Discriminant Analysis (DLDA)
#'
#' Given a set of training data, this function builds the Diagonal Linear
#' Discriminant Analysis (DLDA) classifier, which is often attributed to Dudoit
#' et al. (2002). The DLDA classifier belongs to the family of Naive Bayes
#' classifiers, where the distributions of each class are assumed to be
#' multivariate normal and to share a common covariance matrix.
#'
#' The DLDA classifier is a modification to the well-known LDA classifier, where
#' the off-diagonal elements of the pooled sample covariance matrix are assumed
#' to be zero -- the features are assumed to be uncorrelated. Under multivariate
#' normality, the assumption uncorrelated features is equivalent to the
#' assumption of independent features. The feature-independence assumption is a
#' notable attribute of the Naive Bayes classifier family. The benefit of these
#' classifiers is that they are fast and have much fewer parameters to estimate,
#' especially when the number of features is quite large.
#'
#' The matrix of training observations are given in \code{x}. The rows of \code{x}
#' contain the sample observations, and the columns contain the features for each
#' training observation.
#'
#' The vector of class labels given in \code{y} are coerced to a \code{factor}.
#' The length of \code{y} should match the number of rows in \code{x}.
#'
#' An error is thrown if a given class has less than 2 observations because the
#' variance for each feature within a class cannot be estimated with less than 2
#' observations.
#'
#' The vector, \code{prior}, contains the \emph{a priori} class membership for
#' each class. If \code{prior} is NULL (default), the class membership
#' probabilities are estimated as the sample proportion of observations belonging
#' to each class. Otherwise, \code{prior} should be a vector with the same length
#' as the number of classes in \code{y}. The \code{prior} probabilties should be
#' nonnegative and sum to one.
#'
#' @export
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param prior vector with prior probabilities for each class. If NULL
#' (default), then equal probabilities are used. See details.
#' @return \code{dlda} object that contains the trained DLDA classifier
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' dlda_out <- dlda(Species ~ ., data = iris[train, ])
#' predicted <- predict(dlda_out, iris[-train, -5])$class
#'
#' dlda_out2 <- dlda(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(dlda_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
dlda <- function(x, ...) {
  UseMethod("dlda")
}

#' @rdname dlda
#' @method dlda default
#' @S3method dlda default
dlda.default <- function(x, y, prior = NULL, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- diag_estimates(x = x, y = y, prior = prior, pool = TRUE)

  # Creates an object of type 'dlda' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "dlda"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname dlda
#' @method dlda formula
#' @S3method dlda formula
dlda.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- dlda.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a DLDA classifier object.
#'
#' Summarizes the trained DLDA classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname dlda
#' @method print dlda
#' @S3method print dlda
#' @export
print.dlda <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Sample Size:\n")
  print(x$N)
  cat("Number of Features:\n")
  print(x$p)
  cat("Classes:\n")
  print(x$groups)
  cat("Prior Probabilties:\n")
  print(sapply(x$est, function(z) z$prior))
}

#' DLDA prediction of the class membership of a matrix of new observations.
#'
#' The DLDA classifier is a modification to LDA, where the off-diagonal elements
#' of the pooled sample covariance matrix are set to zero.
#' 
#' @rdname dlda
#' @method predict dlda
#' @S3method predict dlda
#' @export
#'
#' @param object trained DLDA object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#' @param ... additional arguments
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @return list predicted class memberships of each row in newdata
predict.dlda <- function(object, newdata, ...) {
	if (!inherits(object, "dlda"))  {
		stop("object not of class 'dlda'")
	}
	if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }

	scores <- apply(newdata, 1, function(obs) {
		sapply(object$est, function(class_est) {
			with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
		})
	})
	
	if (is.vector(scores)) {
		min_scores <- which.min(scores)
	} else {
		min_scores <- apply(scores, 2, which.min)
	}

	class <- factor(object$groups[min_scores], levels = object$groups)
	
	list(class = class, scores = scores)
}
