#' Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
#'
#' Given a set of training data, this function builds the Shrinkage-based
#' Diagonal Quadratic Discriminant Analysis (SDQDA) classifier, which is based
#' on the DQDA classifier, often attributed to Dudoit et al. (2002). The DQDA
#' classifier belongs to the family of Naive Bayes classifiers, where the
#' distributions of each class are assumed to be multivariate normal. To improve
#' the estimation of the class variances, Pang et al. (2009) proposed the SDQDA
#' classifier which uses a shrinkage-based estimators of each class covariance
#' matrix.
#'
#' The DQDA classifier is a modification to the well-known QDA classifier, where
#' the off-diagonal elements of the pooled covariance matrix are assumed to be
#' zero -- the features are assumed to be uncorrelated. Under multivariate
#' normality, the assumption uncorrelated features is equivalent to the
#' assumption of independent features. The feature-independence assumption is a
#' notable attribute of the Naive Bayes classifier family. The benefit of these
#' classifiers is that they are fast and have much fewer parameters to estimate,
#' especially when the number of features is quite large.
#'
#' The matrix of training observations are given in \code{x}. The rows of
#' \code{x} contain the sample observations, and the columns contain the
#' features for each training observation.
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
#' probabilities are estimated as the sample proportion of observations
#' belonging to each class. Otherwise, \code{prior} should be a vector with the
#' same length as the number of classes in \code{y}. The \code{prior}
#' probabilties should be nonnegative and sum to one.
#'
#' @export
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param prior vector with prior probabilities for each class. If NULL
#' (default), then equal probabilities are used. See details.
#' @param num_alphas the number of values used to find the optimal amount of
#' shrinkage
#' @return \code{sdqda} object that contains the trained SDQDA classifier
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' @examples
#' set.seed(42)
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' sdqda_out <- sdqda(Species ~ ., data = iris[train, ])
#' predicted <- predict(sdqda_out, iris[-train, -5])$class
#'
#' sdqda_out2 <- sdqda(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(sdqda_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
sdqda <- function(x, ...) {
  UseMethod("sdqda")
}

#' @rdname sdqda
#' @method sdqda default
#' @S3method sdqda default
sdqda.default <- function(x, y, prior = NULL, num_alphas = 101, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- diag_estimates(x, y, prior, pool = FALSE)

  # Calculates the shrinkage-based estimator for each diagonal sample class
  # covariance matrix. We add these to the corresponding obj$est$var_shrink
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$var_shrink <- var_shrinkage(
			N = obj$est[[k]]$n,
			K = 1,
			var_feature = obj$est[[k]]$var,
			num_alphas = num_alphas,
			t = -1
		)
  }

  # Creates an object of type 'sdqda' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "sdqda"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname sdqda
#' @method sdqda formula
#' @S3method sdqda formula
sdqda.formula <- function(formula, data, prior = NULL, num_alphas = 101, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- sdqda.default(x = x, y = y, prior = prior, num_alphas = num_alphas)

  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a SDQDA classifier object.
#'
#' Summarizes the trained SDQDA classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname sdqda
#' @method print sdqda
#' @S3method print sdqda
#' @export
print.sdqda <- function(x, ...) {
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

#' SDQDA prediction of the class membership of a matrix of new observations.
#'
#' The SDQDA classifier is a modification to QDA, where the off-diagonal
#' elements of the pooled sample covariance matrix are set to zero. To improve
#' the estimation of the pooled variances, we use a shrinkage method from Pang
#' et al.  (2009).
#' 
#' @rdname sdqda
#' @method predict sdqda
#' @S3method predict sdqda
#' @export
#'
#' @param object trained SDQDA object
#' @param newdata matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @param ... additional arguments
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' @return list predicted class memberships of each row in newdata
predict.sdqda <- function(object, newdata, ...) {
	if (!inherits(object, "sdqda"))  {
		stop("object not of class 'sdqda'")
	}
	if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }

	scores <- apply(newdata, 1, function(obs) {
		sapply(object$est, function(class_est) {
			with(class_est, sum((obs - xbar)^2 / var_shrink + log(var_shrink))
           + log(prior))
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

