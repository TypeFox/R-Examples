#' Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis (SmDQDA) from
#' Tong, Chen, and Zhao (2012)
#'
#' Given a set of training data, this function builds the Shrinkage-mean-based
#' Diagonal Quadratic Discriminant Analysis (SmDQDA) classifier from Tong, Chen,
#' and Zhao (2012). The SmDQDA classifier incorporates a Lindley-type shrunken
#' mean estimator into the DQDA classifier from Dudoit et al. (2002). For more
#' about the DQDA classifier, see \code{\link{dqda}}.
#'
#' The DQDA classifier is a modification to the well-known QDA classifier, where
#' the off-diagonal elements of each class covariance matrix are assumed
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
#' @return \code{smdqda} object that contains the trained SmDQDA classifier
#'
#' @references Tong, T., Chen, L., and Zhao, H. (2012), "Improved Mean
#' Estimation and Its Application to Diagonal Discriminant Analysis,"
#' Bioinformatics, 28, 4, 531-537.
#' \url{http://bioinformatics.oxfordjournals.org/content/28/4/531.long}
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' smdqda_out <- smdqda(Species ~ ., data = iris[train, ])
#' predicted <- predict(smdqda_out, iris[-train, -5])$class
#'
#' smdqda_out2 <- smdqda(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(smdqda_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
smdqda <- function(x, ...) {
  UseMethod("smdqda")
}

#' @rdname smdqda
#' @method smdqda default
#' @S3method smdqda default
smdqda.default <- function(x, y, prior = NULL, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  obj <- diag_estimates(x = x, y = y, prior = prior, pool = FALSE,
                        est_mean = "tong")

  # Creates an object of type 'smdqda' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "smdqda"
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname smdqda
#' @method smdqda formula
#' @S3method smdqda formula
smdqda.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- smdqda.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a SmDQDA classifier object.
#'
#' Summarizes the trained SmDQDA classifier in a nice manner.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname smdqda
#' @method print smdqda
#' @S3method print smdqda
#' @export
print.smdqda <- function(x, ...) {
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

#' SmDQDA prediction of the class membership of a matrix of new observations.
#'
#' The SmDQDA classifier is a modification to QDA, where the off-diagonal elements
#' of the pooled sample covariance matrix are set to zero.
#' 
#' @rdname smdqda
#' @method predict smdqda
#' @S3method predict smdqda
#' @export
#'
#' @param object trained SmDQDA object
#' @param newdata matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @param ... additional arguments
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @return list predicted class memberships of each row in newdata
predict.smdqda <- function(object, newdata, ...) {
	if (!inherits(object, "smdqda"))  {
		stop("object not of class 'smdqda'")
	}
	if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }

	scores <- apply(newdata, 1, function(obs) {
		sapply(object$est, function(class_est) {
			with(class_est, sum((obs - xbar)^2 / var + log(var)) + log(prior))
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
