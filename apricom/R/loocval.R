#' Leave-one-out Cross-validation-derived Shrinkage
#'
#' Shrink regression coefficients using a shrinkage factor derived using
#' leave-one-out cross-validation.
#'
#' This function applies leave-one-out cross-validation to a dataset in order to
#' derive a shrinkage factor and apply it to the regression coefficients. One row
#' of the data is used as a validation set, while the remaining data is used as
#' a training set. Regression coefficients are estimated in the training set, and
#' then a shrinkage factor is estimated using the validation set. This process is
#' repeated so that each data row is used as the validation set once. The mean of
#' the shrinkage factors is then applied to the original regression coeffients,
#' and the regression intercept may be re-estimated. This process may be repeated
#' nreps times but each rep should yield the same shrunken coefficients.
#'
#' This process can currently be applied to linear or logistic regression models.
#'
#' @importFrom stats glm.fit binomial
#'
#' @param dataset   a dataset for regression analysis. Data should be in the form
#' of a matrix, with the outcome variable as the final column. Application of the
#' \code{\link{datashape}} function beforehand is recommended, especially if
#'  categorical predictors are present. For regression with an intercept included
#'  a column vector of 1s should be included before the dataset (see examples)
#' @param model    type of regression model. Either "linear" or "logistic".
#' @param nreps    the number of times to replicate the cross-validation process.
#' @param sdm      a shrinkage design matrix.
#' @param int      logical. If TRUE the model will include a
#'                 regression intercept.
#' @param int.adj  logical. If TRUE the regression intercept will be
#'                 re-estimated after shrinkage of the regression coefficients.
#'
#' @return \code{loocval} returns a list containing the following:
#' @return \item{raw.coeff}{the raw regression model coefficients,
#'               pre-shrinkage.}
#' @return \item{shrunk.coeff}{the shrunken regression model coefficients.}
#' @return \item{lambda}{the mean shrinkage factor over nreps
#'               cross-validation replicates.}
#' @return \item{nreps}{the number of cross-validation replicates.}
#' @return \item{sdm}{the shrinkage design matrix used to apply
#'              the shrinkage factor(s) to the regression coefficients.}
#'
#' @note \strong{Warning:} this method is not recommended for use in practice. Due to the
#'                         high variance and inherent instability of leave-one-out methods
#'                         the value of the shrinkage factor may be extreme.
#'
#' @examples
#'## Example 1: Linear regression using the iris dataset
#'## Leave-one-out cross-validation-derived shrinkage
#' data(iris)
#' iris.data <- as.matrix(iris[, 1:4])
#' iris.data <- cbind(1, iris.data)
#' sdm1 <- matrix(c(0, 1, 1, 1), nrow = 1)
#' set.seed(123)
#' loocval(dataset = iris.data, model = "linear", sdm = sdm1,
#' int = TRUE, int.adj = TRUE)
#'
#'## Example 2: logistic regression using a subset of the mtcars data
#'## Leave-one-out cross-validation-derived shrinkage
#' data(mtcars)
#' mtc.data <- cbind(1,datashape(mtcars, y = 8, x = c(1, 6, 9)))
#' head(mtc.data)
#' set.seed(123)
#' loocval(dataset = mtc.data, model = "logistic")

loocval<- function(dataset, model, nreps = 1, sdm, int = TRUE, int.adj) {

  nc <- dim(dataset)[2]
  nr <- dim(dataset)[1]
  lambda <- c(rep(0, nreps))
  B.shrunk <- matrix(c(rep(0, nreps * (nc - 1))), ncol = nreps)

  if (model == "linear") B <- ols.rgr(dataset)
  else B <- ml.rgr(dataset)

  if (int == FALSE) int.adj <- FALSE
  if (int == TRUE & missing(int.adj)) int.adj <- TRUE
    if (missing(sdm)){
      if (int == "FALSE") sdm <- matrix(rep(1, nc - 1), nrow = 1)
      else sdm <- matrix(c(0, rep(1, nc - 2)), nrow = 1)
    }

  sdm2 <- t(sdm)

  for (i in 1:nreps) {
    s <- c(rep(0, nr))
    ### LINEAR REGRESSION MODEL
    if (model == "linear") {
      for (k in 1:nr) {
        D.val <- matrix(dataset[k, ], nrow = 1)
        D.train <- dataset[-k, ]
        train.model <- ols.rgr(D.train)
        # Estimate shrinkage factor
        s[k] <- c(ols.shrink(train.model, D.val, sdm))
      }
    lambda[i] <- mean(s)
      # Apply shrinkage factor
      B.shrunk[, i] <- matrix(diag(as.vector(B)) %*% sdm2 %*% lambda[i] +
                  diag(as.vector(B)) %*% apply(1 - sdm2, 1, min), ncol = 1)

      if (int.adj == TRUE) {
        # re-estimate the intercept
        new.int <- mean((dataset[, nc]) - ((dataset[, 2:(nc - 1)]) %*%
                        B.shrunk[-1, i]))
        B.shrunk[, i] <- c(new.int, B.shrunk[-1, i])
      }


      ### LOGISTIC REGRESSION MODEL
    } else {

      for (k in 1:nr) {
        D.val <- matrix(dataset[k, ], nrow = 1)
        D.train <- dataset[-k, ]
        train.model <- ols.rgr(D.train)
        # Estimate shrinkage factor
        s[k] <- c(ml.shrink(train.model, D.val))
      }
      lambda[i] <- mean(s)
      # Apply shrinkage factor
      B.shrunk[, i] <- matrix(diag(as.vector(B)) %*% sdm2 %*% lambda[i] +
                  diag(as.vector(B)) %*% apply(1 - sdm2, 1, min), ncol = 1)

      if (int.adj == TRUE) {
        # re-estimate the intercept
        dataset <- as.matrix(dataset)
        X.i <-  matrix(dataset[, 1], ncol = 1)
        Y <- dataset[, nc]
        offs <- as.vector((dataset[, 2:(nc - 1)]) %*% B.shrunk[-1, i])
        new.int <-glm.fit(X.i, Y, family = binomial(link = "logit"),
                          offset = offs)$coefficients

        B.shrunk[, i] <- c(new.int, B.shrunk[-1, i])
      }
    }
  }
  if (nreps > 1) {
  B.shrunk.out <- rowMeans(B.shrunk)
  lambda.out <- mean(lambda)
  return(list(raw.coeff = B, shrunk.coeff = matrix(B.shrunk.out, ncol = 1),
              lambda = lambda.out, nRounds = nreps, sdm = sdm))
  } else {
    return(list(raw.coeff = B, shrunk.coeff = B.shrunk, lambda = lambda,
                nRounds = nreps, sdm = sdm))
  }
}

