#' The Minimum Distance Rule using Moore-Penrose Inverse (MDMP) classifier
#'
#' Given a set of training data, this function builds the MDMP classifier from
#' Srivistava and Kubokawa (2007). The MDMP classifier is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Srivastava and Kubokawa (2007) have
#' proposed a modification of the standard maximum likelihood estimator of the
#' pooled covariance matrix, where only the largest 95% of the eigenvalues and
#' their corresponding eigenvectors are kept. The value of 95% is the default
#' and can be changed via the \code{eigen_pct} argument.
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
#' @param eigen_pct the percentage of eigenvalues kept
#' @return \code{mdmp} object that contains the trained MDMP classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' mdmp_out <- mdmp(Species ~ ., data = iris[train, ])
#' predicted <- predict(mdmp_out, iris[-train, -5])$class
#'
#' mdmp_out2 <- mdmp(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(mdmp_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
mdmp <- function(x, ...) {
  UseMethod("mdmp")
}

#' @rdname mdmp
#' @method mdmp default
#' @S3method mdmp default
mdmp.default <- function(x, y, prior = NULL, eigen_pct = 0.95, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)

  # Removes original pooled covariance matrix to reduce memory usage 
  obj$cov_pool <- NULL

  # trace(cov_kept) / trace(cov_pool) \approx eigen_pct
  # as described in the middle of page 125
  kept_evals <- with(cov_eigen,
                     which(cumsum(values) / sum(values) < eigen_pct))

  # Computes the pseudoinverse of the resulting covariance matrix estimator
  evals_inv <- 1 / cov_eigen$values[kept_evals]
  obj$cov_inv <- with(cov_eigen,
                      tcrossprod(vectors[, kept_evals] %*% diag(evals_inv),
                                 vectors[, kept_evals]))

  # Creates an object of type 'mdmp' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "mdmp"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname mdmp
#' @method mdmp formula
#' @S3method mdmp formula
mdmp.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- mdmp.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a MDMP classifier object.
#'
#' Summarizes the trained mdmp classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname mdmp
#' @method print mdmp
#' @S3method print mdmp
#' @export
print.mdmp <- function(x, ...) {
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
#' Minimum Distance Rule using Moore-Penrose Inverse (MDMP) classifier
#'
#' The MDMP classifier from Srivistava and Kubokawa (2007) is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Srivastava and Kubokawa (2007) have
#' proposed a modification of the standard maximum likelihood estimator of the
#' pooled covariance matrix, where only the largest 95% of the eigenvalues and
#' their corresponding eigenvectors are kept.
#' 
#' @rdname mdmp
#' @method predict mdmp
#' @S3method predict mdmp
#' @export
#'
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
#' @param object trained mdmp object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#' @param ... additional arguments
#' @return list predicted class memberships of each row in newdata
predict.mdmp <- function(object, newdata, ...) {
	if (!inherits(object, "mdmp"))  {
		stop("object not of class 'mdmp'")
	}
	if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }

  # Calculates the discriminant scores for each test observation
	scores <- apply(newdata, 1, function(obs) {
		sapply(object$est, function(class_est) {
			with(class_est, quadform(object$cov_inv, obs - xbar) + log(prior))
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
