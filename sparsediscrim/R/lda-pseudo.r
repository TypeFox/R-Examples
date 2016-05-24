#' Linear Discriminant Analysis (LDA) with the Moore-Penrose Pseudo-Inverse
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix.
#' When the pooled sample covariance matrix is singular, the linear discriminant
#' function is incalculable. A common method to overcome this issue is to
#' replace the inverse of the pooled sample covariance matrix with the
#' Moore-Penrose pseudo-inverse, which is unique and always exists. Note that
#' when the pooled sample covariance matrix is nonsingular, it is equal to the
#' pseudo-inverse.
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
#' @param tol tolerance value below which eigenvalues are considered numerically
#' equal to 0
#' @return \code{lda_pseudo} object that contains the trained lda_pseudo
#' classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' lda_pseudo_out <- lda_pseudo(Species ~ ., data = iris[train, ])
#' predicted <- predict(lda_pseudo_out, iris[-train, -5])$class
#'
#' lda_pseudo_out2 <- lda_pseudo(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(lda_pseudo_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
lda_pseudo <- function(x, ...) {
  UseMethod("lda_pseudo")
}

#' @rdname lda_pseudo
#' @method lda_pseudo default
#' @S3method lda_pseudo default
lda_pseudo.default <- function(x, y, prior = NULL, tol = 1e-8, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Calculates the Moore-Penrose pseudo inverse based on the pooled sample
  # covariance matrix
  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)
  evals <- cov_eigen$values
  evals_inv <- rep.int(0, times = length(evals))
  evals_inv[evals > tol] <- 1 / evals[evals > tol]

  # Removes original pooled covariance matrix to reduce memory usage 
  obj$cov_pool <- NULL

  obj$cov_inv <- with(cov_eigen,
                      tcrossprod(vectors %*% diag(evals_inv), vectors))

  # Creates an object of type 'lda_pseudo' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "lda_pseudo"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname lda_pseudo
#' @method lda_pseudo formula
#' @S3method lda_pseudo formula
lda_pseudo.formula <- function(formula, data, prior = NULL, tol = 1e-8, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_pseudo.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a lda_pseudo classifier object.
#'
#' Summarizes the trained lda_pseudo classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname lda_pseudo
#' @method print lda_pseudo
#' @S3method print lda_pseudo
#' @export
print.lda_pseudo <- function(x, ...) {
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
#' Discriminant Analysis (LDA) with the Moore-Penrose Pseudo-Inverse
#'
#' The Linear Discriminant Analysis (LDA) classifier involves the assumption
#' that the distributions of each class are assumed to be multivariate normal
#' and share a common covariance matrix. When the pooled sample covariance
#' matrix is singular, the linear discriminant function is incalculable. A
#' common method to overcome this issue is to replace the inverse of the pooled
#' sample covariance matrix with the Moore-Penrose pseudo-inverse, which is
#' unique and always exists. Note that when the pooled sample covariance matrix
#' is nonsingular, it is equal to the pseudo-inverse.
#' 
#' @rdname lda_pseudo
#' @method predict lda_pseudo
#' @S3method predict lda_pseudo
#'
#' @export
#'
#' @param object trained lda_pseudo object
#' @param newdata matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @param ... additional arguments
#' @return list predicted class memberships of each row in newdata
predict.lda_pseudo <- function(object, newdata, ...) {
	if (!inherits(object, "lda_pseudo"))  {
		stop("object not of class 'lda_pseudo'")
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

