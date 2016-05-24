#' The Minimum Distance Empirical Bayesian Estimator (MDEB) classifier
#'
#' Given a set of training data, this function builds the MDEB classifier from
#' Srivistava and Kubokawa (2007). The MDEB classifier is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Rather than using the standard maximum
#' likelihood estimator of the pooled covariance matrix, Srivastava and Kubokawa
#' (2007) have proposed an Empirical Bayes estimator where the eigenvalues of
#' the pooled sample covariance matrix are shrunken towards the identity matrix:
#' the shrinkage constant has a closed form and is quick to calculate.
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
#' @return \code{mdeb} object that contains the trained MDEB classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' mdeb_out <- mdeb(Species ~ ., data = iris[train, ])
#' predicted <- predict(mdeb_out, iris[-train, -5])$class
#'
#' mdeb_out2 <- mdeb(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(mdeb_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
mdeb <- function(x, ...) {
  UseMethod("mdeb")
}

#' @rdname mdeb
#' @method mdeb default
#' @S3method mdeb default
mdeb.default <- function(x, y, prior = NULL, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Creates an object of type 'mdeb' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "mdeb"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname mdeb
#' @method mdeb formula
#' @S3method mdeb formula
mdeb.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- mdeb.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a MDEB classifier object.
#'
#' Summarizes the trained mdeb classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname mdeb
#' @method print mdeb
#' @S3method print mdeb
#' @export
print.mdeb <- function(x, ...) {
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

#' Predicts of class membership of a matrix of new observations using the
#' Minimum Distance Empirical Bayesian Estimator (MDEB) classifier
#'
#' The MDEB classifier from Srivistava and Kubokawa (2007) is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Rather than using the standard maximum
#' likelihood estimator of the pooled covariance matrix, Srivastava and Kubokawa
#' (2007) have proposed an Empirical Bayes estimator where the eigenvalues of
#' the pooled sample covariance matrix are shrunken towards the identity matrix:
#' the shrinkage constant has a closed form and is quick to calculate
#' 
#' @rdname mdeb
#' @method predict mdeb
#' @S3method predict mdeb
#' @export
#'
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
#' @param object trained mdeb object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#' @param ... additional arguments
#' @return list predicted class memberships of each row in newdata
predict.mdeb <- function(object, newdata, ...) {
	if (!inherits(object, "mdeb"))  {
		stop("object not of class 'mdeb'")
	}
	if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }

  # Calculates the MDEB shrinkage constant and then computes the inverse of the
  # MDEB covariance matrix estimator
  shrink <- with(object, sum(diag(cov_pool)) / min(N, p))
  cov_inv <- with(object, solve(cov_pool + shrink * diag(p)))

  # Calculates the discriminant scores for each test observation
	scores <- apply(newdata, 1, function(obs) {
		sapply(object$est, function(class_est) {
			with(class_est, quadform(cov_inv, obs - xbar) + log(prior))
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

