#' Cross-validation-derived Shrinkage After Estimation
#'
#' Shrink regression coefficients using a Cross-validation-derived
#' shrinkage factor.
#'
#' This function applies k-fold cross-validation to a dataset in order to
#' derive a shrinkage factor and apply it to the regression coefficients.
#' Data is randomly partitioned into k equally sized sets. One set is used as a
#' validation set, while the remaining data is used as a training set. Regression
#' coefficients are estimated in the training set, and then a shrinkage factor
#' is estimated using the validation set. This process is repeated so that each
#' partitioned set is used as the validation set once, resulting in k folds.
#' The mean of k shrinkage factors is then applied to the original
#' regression coeffients, and the regression intercept may be re-estimated.
#' This process is repeated nreps times and the mean of the regression
#' coefficients is returned.
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
#' @param k        the number of cross-validation folds. This number must be within
#'                 the range 1 < k <= 0.5 * number of observations
#' @param nreps    the number of times to replicate the cross-validation process.
#' @param sdm      a shrinkage design matrix. For examples, see \code{\link{ols.shrink}}
#' @param int      logical. If TRUE the model will include a
#'                 regression intercept.
#' @param int.adj  logical. If TRUE the regression intercept will be
#'                 re-estimated after shrinkage of the regression coefficients.
#'
#' @return \code{kcrossval} returns a list containing the following:
#' @return \item{raw.coeff}{the raw regression model coefficients,
#'               pre-shrinkage.}
#' @return \item{shrunk.coeff}{the shrunken regression model coefficients.}
#' @return \item{lambda}{the mean shrinkage factor over nreps
#'               cross-validation replicates.}
#' @return \item{nFolds}{the number of cross-validation folds.}
#' @return \item{nreps}{the number of cross-validation replicates.}
#' @return \item{sdm}{the shrinkage design matrix used to apply
#'              the shrinkage factor(s) to the regression coefficients.}
#'
#' @examples
#'## Example 1: Linear regression using the iris dataset
#'## 2-fold Cross-validation-derived shrinkage with 50 reps
#' data(iris)
#' iris.data <- as.matrix(iris[, 1:4])
#' iris.data <- cbind(1, iris.data)
#' sdm1 <- matrix(c(0, 1, 1, 1), nrow = 1)
#' kcrossval(dataset = iris.data, model = "linear", k = 2,
#' nreps = 50, sdm = sdm1, int = TRUE, int.adj = TRUE)
#'
#'## Example 2: logistic regression using a subset of the mtcars data
#'## 10-fold CV-derived shrinkage (uniform shrinkage and intercept re-estimation)
#' data(mtcars)
#' mtc.data <- cbind(1,datashape(mtcars, y = 8, x = c(1, 6, 9)))
#' head(mtc.data)
#' set.seed(321)
#' kcrossval(dataset = mtc.data, model = "logistic",
#' k = 10, nreps = 10)


kcrossval<- function(dataset, model, k, nreps, sdm, int = TRUE, int.adj) {

  nc <- dim(dataset)[2]
  lambda <- c(rep(0, nreps))
  B.shrunk <- matrix(c(rep(0, nreps * (nc - 1))), ncol = nreps)

  if (model == "linear") B <- ols.rgr(dataset)
  else B <- ml.rgr(dataset)

  if (int == FALSE) int.adj <- FALSE
  if (int == TRUE & missing(int.adj)) int.adj <- TRUE
  if (missing(nreps)) nreps <- 1
  if (missing(k)) stop("Number of folds (k) for cross-validation must be specified")
  if (missing(sdm)){
    if (int == "FALSE") sdm <- matrix(rep(1, nc - 1), nrow = 1)
    else sdm <- matrix(c(0, rep(1, nc - 2)), nrow = 1)
  }

  sdm2 <- t(sdm)

  for (i in 1:nreps) {

    Dats <- grandpart(dataset, k)
    s <- c(rep(0, k))

### LINEAR REGRESSION MODEL
    if (model == "linear") {

      for (l in 1:k) {
        if (k == 2) D.train <- as.matrix(Dats[[-l]])
        else D.train <- do.call("rbind", Dats[-l])
        D.val <- as.matrix(Dats[[l]])
        train.model <- ols.rgr(D.train)
# Estimate shrinkage factor
        s[l] <- c(ols.shrink(train.model, D.val, sdm))
      }

      lambda[i] <- mean(s)
# Apply shrinkage factor
      B.shrunk[, i] <- matrix(diag(as.vector(B)) %*% sdm2 %*% lambda[i] + diag(as.vector(B)) %*% apply(1 - sdm2, 1, min), ncol = 1)

      if (int.adj == TRUE) {
# re-estimate the intercept
        new.int <- mean((dataset[, nc]) - ((dataset[, 2:(nc - 1)]) %*% B.shrunk[-1, i]))
        B.shrunk[, i] <- c(new.int, B.shrunk[-1, i])
      }


### LOGISTIC REGRESSION MODEL
    } else {

        for (l in 1:k) {
          if (k == 2) D.train <- as.matrix(Dats[[-l]])
          else D.train <- do.call("rbind", Dats[-l])
          D.val <- as.matrix(Dats[[l]])
          train.model <- ml.rgr(D.train)
# Estimate shrinkage factor
          s[l] <- c(ml.shrink(train.model, D.val))
        }

        lambda[i] <- mean(s)
# Apply shrinkage factor
        B.shrunk[, i] <- matrix(diag(as.vector(B)) %*% sdm2 %*% lambda[i] + diag(as.vector(B)) %*% apply(1 - sdm2, 1, min), ncol = 1)

        if (int.adj == TRUE) {
# re-estimate the intercept
          dataset <- as.matrix(dataset)
          X.i <-  matrix(dataset[, 1], ncol = 1)
          Y <- dataset[, nc]
          offs <- as.vector((dataset[, 2:(nc - 1)]) %*% B.shrunk[-1, i])
          new.int <-glm.fit(X.i, Y, family = binomial(link = "logit"), offset = offs)$coefficients

          B.shrunk[, i] <- c(new.int, B.shrunk[-1, i])
        }
    }
  }

    B.shrunk.out <- rowMeans(B.shrunk)
    lambda.out <- mean(lambda)
    return(list(raw.coeff = B, shrunk.coeff = matrix(B.shrunk.out, ncol = 1), lambda = lambda.out, nFolds = k, nRounds = nreps, sdm = sdm))

}

