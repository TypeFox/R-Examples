#' Bootstrap-derived Shrinkage After Estimation
#'
#' Shrink regression coefficients using a bootstrap-derived shrinkage factor.
#'
#' This function applies bootstrapping to a dataset in order to derive a shrinkage
#' factor and apply it to the regression coefficients. Regression coefficients are
#' estimated in a bootstrap sample, and then a shrinkage factor is estimated using
#' the input data. The mean of N shrinkage factors is then applied to the original
#' regression coeffients, and the regression intercept may be re-estimated.
#'
#' This process can currently be applied to linear or logistic regression models.
#'
#' @importFrom stats glm.fit binomial
#'
#' @param dataset   a dataset for regression analysis. Data should be in the form
#'                  of a matrix, with the outcome variable as the final column.
#'                  Application of the \code{\link{datashape}} function beforehand
#'                  is recommended, especially if categorical predictors are present.
#'                  For regression with an intercept included a column vector of
#'                  1s should be included before the dataset (see examples).
#' @param model     type of regression model. Either "linear" or "logistic".
#' @param N         the number of times to replicate the bootstrapping process
#' @param sdm       a shrinkage design matrix. For examples, see \code{\link{ols.shrink}}
#' @param int       logical. If TRUE the model will include a regression intercept.
#' @param int.adj   logical. If TRUE the regression intercept will be re-estimated
#'                  after shrinkage of the regression coefficients.
#'
#' @return \code{bootval} returns a list containing the following:
#' @return \item{raw.coeff}{the raw regression model coefficients, pre-shrinkage.}
#' @return \item{shrunk.coeff}{the shrunken regression model coefficients}
#' @return \item{lambda}{the mean shrinkage factor over N bootstrap replicates}
#' @return \item{N}{the number of bootstrap replicates}
#' @return \item{sdm}{the shrinkage design matrix used to apply the shrinkage
#'   factor(s) to the regression coefficients}
#'
#' @examples
#'## Example 1: Linear regression using the iris dataset
#' data(iris)
#' iris.data <- as.matrix(iris[, 1:4])
#' iris.data <- cbind(1, iris.data)
#' sdm1 <- matrix(c(0, 1, 1, 1), nrow = 1)
#' set.seed(777)
#' bootval(dataset = iris.data, model = "linear", N = 200, sdm = sdm1,
#' int = TRUE, int.adj = TRUE)
#'
#'## Example 2: logistic regression using a subset of the mtcars data
#' data(mtcars)
#' mtc.data <- cbind(1,datashape(mtcars, y = 8, x = c(1, 6, 9)))
#' head(mtc.data)
#' set.seed(777)
#' bootval(dataset = mtc.data, model = "logistic", N = 500)

bootval <- function(dataset, model, N, sdm, int = TRUE, int.adj) {

  m <- dim(dataset)[1]
  nc <- dim(dataset)[2]
  s <- c(rep(0, N))

  if (missing(sdm)){
    if (int == FALSE) sdm <- matrix(rep(1, nc - 1), nrow = 1)
    else sdm <- matrix(c(0, rep(1, nc - 2)), nrow = 1)
  }

  if (int == FALSE) int.adj <- FALSE
  if (int == TRUE & missing(int.adj)) int.adj <- TRUE


### LINEAR REGRESSION MODEL

  if (model == "linear") {

# Bootstrap estimation of a shrinkage factor
    for(k in 1:N) {
      indces <- sample(m, replace = TRUE)
      dboot <- dataset[indces, ]
      b <- ols.rgr(dboot)
      s[k] <- ols.shrink(b, dataset, sdm)
    }

    lambda <- mean(s)

# Apply shrinkage factor
    betas <- ols.rgr(dataset)
    sdm <- t(sdm)
    B.shrunk <- matrix(diag(as.vector(betas)) %*% sdm %*% lambda +
                       diag(as.vector(betas)) %*% apply(1 - sdm, 1, min), ncol = 1)

  if (int.adj == TRUE) {

# re-estimate the intercept
      new.int <- mean((dataset[, nc]) - ((dataset[, 2:(nc - 1)]) %*% B.shrunk[-1]))

      B.shrunk <- c(new.int, B.shrunk[-1])
    }

    return(list(raw.coeff = betas, shrunk.coeff = matrix(B.shrunk, ncol = 1),
                lambda = lambda, N = N, sdm = t(sdm)))

### LOGISTIC REGRESSION MODEL

    } else {

# Bootstrap estimation of a shrinkage factor
      for(k in 1:N) {
        indces <- sample(m, replace = TRUE)
        dboot <- dataset[indces, ]
        b <- ml.rgr(dboot)
        s[k] <- ml.shrink(b, dataset)  ## CONSIDER INC. METHOD FOR OPTIM

      }

      lambda <- mean(s)

# Apply shrinkage factor
      betas <- ml.rgr(dataset)
      sdm <- t(sdm)
      B.shrunk <- matrix(diag(as.vector(betas)) %*% sdm %*% lambda +
                      diag(as.vector(betas)) %*% apply(1 - sdm, 1, min), ncol = 1)

      if (int.adj == TRUE) {

# re-estimate the intercept

        dataset <- as.matrix(dataset)
        X.i <-  matrix(dataset[, 1], ncol = 1)
        Y <- dataset[, nc]
        offs <- as.vector((dataset[, 2:(nc - 1)]) %*% B.shrunk[-1])
        new.int <-glm.fit(X.i, Y, family = binomial(link = "logit"),
                          offset = offs)$coefficients

        B.shrunk <- c(new.int, B.shrunk[-1])
      }

      return(list(raw.coeff = betas, shrunk.coeff = matrix(B.shrunk, ncol = 1),
                  lambda = lambda, N = N, sdm = t(sdm)))
    }
}
