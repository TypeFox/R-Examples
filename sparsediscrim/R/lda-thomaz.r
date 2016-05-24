#' Linear Discriminant Analysis using the Thomaz-Kitani-Gillies Covariance
#' Matrix Estimator
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix. When the
#' pooled sample covariance matrix is singular, the linear discriminant function
#' is incalculable. This function replaces the pooled sample covariance matrix
#' with a regularized estimator from Thomaz et al. (2006), where the smallest
#' eigenvalues are replaced with the average eigenvalue. Specifically, small
#' eigenvalues here means that the eigenvalues are less than the average
#' eigenvalue.
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
#' @return \code{lda_thomaz} object that contains the trained classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' lda_thomaz_out <- lda_thomaz(Species ~ ., data = iris[train, ])
#' predicted <- predict(lda_thomaz_out, iris[-train, -5])$class
#'
#' lda_thomaz_out2 <- lda_thomaz(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(lda_thomaz_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
#' @references Thomaz, C. E., Kitani, E. C., and Gillies, D. F. (2006). "A
#' maximum uncertainty LDA-based approach for limited sample size problems with
#' application to face recognition," J. Braz. Comp. Soc., 12, 2, 7-18.
lda_thomaz <- function(x, ...) {
  UseMethod("lda_thomaz")
}

#' @rdname lda_thomaz
#' @method lda_thomaz default
#' @S3method lda_thomaz default
lda_thomaz.default <- function(x, y, prior = NULL, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Computes eigenvalue decomposition of pooled sample covariance matrix
  # Then regularizes the estimator based on Thomaz et al.'s (2006) method
  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)
  evals <- cov_eigen$values
  mean_eval <- mean(evals)
  evals[evals < mean_eval] <- mean_eval

  # Removes original pooled covariance matrix to reduce memory usage 
  obj$cov_pool <- NULL

  # Computes the inverse of the resulting covariance matrix estimator
  obj$cov_inv <- with(cov_eigen,
                      tcrossprod(vectors %*% diag(1 / evals), vectors))

  # Creates an object of type 'lda_thomaz' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "lda_thomaz"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname lda_thomaz
#' @method lda_thomaz formula
#' @S3method lda_thomaz formula
lda_thomaz.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_thomaz.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a lda_thomaz classifier object.
#'
#' Summarizes the trained lda_thomaz classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname lda_thomaz
#' @method print lda_thomaz
#' @S3method print lda_thomaz
#' @export
print.lda_thomaz <- function(x, ...) {
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

#' Predicts of class membership of a matrix of new observations using Linear
#' Discriminant Analysis (LDA) using the Schafer-Strimmer Covariance Matrix
#' Estimator
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix. When the
#' pooled sample covariance matrix is singular, the linear discriminant function
#' is incalculable. This function replaces the pooled sample covariance matrix
#' with a regularized estimator from Thomaz et al. (2006), where the smallest
#' eigenvalues are replaced with the average eigenvalue. Specifically, small
#' eigenvalues here means that the eigenvalues are less than the average
#' eigenvalue.
#' 
#' @rdname lda_thomaz
#' @method predict lda_thomaz
#' @S3method predict lda_thomaz
#' @export
#'
#' @references Thomaz, C. E., Kitani, E. C., and Gillies, D. F. (2006). "A
#' maximum uncertainty LDA-based approach for limited sample size problems with
#' application to face recognition," J. Braz. Comp. Soc., 12, 2, 7-18.
#' @param object trained lda_thomaz object
#' @param newdata matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @param ... additional arguments
#' @return list predicted class memberships of each row in newdata
predict.lda_thomaz <- function(object, newdata, ...) {
	if (!inherits(object, "lda_thomaz"))  {
		stop("object not of class 'lda_thomaz'")
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

