#' Split-sample-derived Shrinkage After Estimation
#'
#' Shrink regression coefficients using a split-sample-derived shrinkage factor.
#'
#' This function applies sample-splitting to a dataset in order to derive a
#' shrinkage factor and apply it to the regression coefficients. Data are
#' randomly split into two sets, a training set and a test set. Regression
#' coefficients are estimated using the training sample, and then a shrinkage
#' factor is estimated using the test set. The mean of N shrinkage factors is
#' then applied to the original regression coeffients, and the regression
#' intercept may be re-estimated.
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
#' @param model     type of regression model. Either "linear" or "logistic".
#' @param nrounds   the number of times to replicate the sample splitting process.
#' @param fract     the fraction of observations designated to the training set
#' @param sdm       a shrinkage design matrix. For examples, see \code{\link{ols.shrink}}
#' @param int       logical. If TRUE the model will include a regression intercept.
#' @param int.adj   logical. If TRUE the regression intercept will be re-estimated
#'                 after shrinkage of the regression coefficients.
#'
#' @return \code{splitval} returns a list containing the following:
#' @return \item{raw.coeff}{the raw regression model coefficients, pre-shrinkage.}
#' @return \item{shrunk.coeff}{the shrunken regression model coefficients}
#' @return \item{lambda}{the mean shrinkage factor over Nrounds split-sample
#'               replicates}
#' @return \item{Nrounds}{the number of rounds of sample splitting}
#' @return \item{sdm}{the shrinkage design matrix used to apply the shrinkage
#'               factor(s) to the regression coefficients}
#'
#' @examples
#'## Example 1: Linear regression using the iris dataset
#'## Split-sample-derived shrinkage with 100 rounds of sample-splitting
#' data(iris)
#' iris.data <- as.matrix(iris[, 1:4])
#' iris.data <- cbind(1, iris.data)
#' sdm1 <- matrix(c(0, 1, 1, 1), nrow = 1)
#' set.seed(321)
#' splitval(dataset = iris.data, model = "linear", nrounds = 100,
#' fract = 0.75, sdm = sdm1, int = TRUE, int.adj = TRUE)
#'
#'## Example 2: logistic regression using a subset of the mtcars data
#'## Split-sample-derived shrinkage
#' data(mtcars)
#' mtc.data <- cbind(1,datashape(mtcars, y = 8, x = c(1, 6, 9)))
#' head(mtc.data)
#' set.seed(123)
#' splitval(dataset = mtc.data, model = "logistic",
#' nrounds = 100, fract = 0.5)

splitval <- function(dataset, model, nrounds, fract, sdm, int = TRUE, int.adj) {

  nr <- dim(dataset)[1]
  nc <- dim(dataset)[2]
  s <- c(rep(0, nrounds))
  if (model == "linear") B <- ols.rgr(dataset)
  else B <- ml.rgr(dataset)

  if (missing(fract)) {fract <- 0.5}
  if (missing(nrounds)) {nrounds <- 1}
  if (missing(sdm)){
    if (int == "FALSE") sdm <- matrix(rep(1, nc - 1), nrow = 1)
    else sdm <- matrix(c(0, rep(1, nc - 2)), nrow = 1)
  }
  if ((nr - (fract * nr)) <= 1) stop("Partitioning fraction is too large")
  if (nr * fract < (nc - 1)) stop("Partitioning fraction is too small")
  if (int == FALSE) int.adj <- FALSE
  if (int == TRUE & missing(int.adj)) int.adj <- TRUE


  for (i in 1:nrounds) {

    if (model == "linear") {
      Dats <- randpart(dataset, fract)
      bf <- ols.rgr(Dats[[1]])
      s[i] <- ols.shrink(bf, Dats[[2]], sdm)

      } else {
        Dats <- randpart(dataset, fract)
        bf <- ml.rgr(Dats[[1]])
        s[i] <- ml.shrink(bf, Dats[[2]])
      }
  }

  lambda <- mean(s)

  if (model == "linear") {

    # Apply shrinkage factor
    B <- ols.rgr(dataset)
    sdm <- t(sdm)
    B.shrunk <- matrix(diag(as.vector(B)) %*% sdm %*% lambda +
                         diag(as.vector(B)) %*% apply(1 - sdm, 1, min), ncol = 1)

    if (int.adj == TRUE) {

      # re-estimate the intercept
      new.int <- mean((dataset[, nc]) - ((dataset[, 2:(nc - 1)]) %*% B.shrunk[-1]))

      B.shrunk <- c(new.int, B.shrunk[-1])
    }

    return(list(raw.coeff = B, shrunk.coeff = matrix(B.shrunk, ncol = 1),
                lambda = lambda, Nrounds = nrounds, sdm = t(sdm)))

  } else {

    # Apply shrinkage factor
    B <- ml.rgr(dataset)
    sdm <- t(sdm)
    B.shrunk <- matrix(diag(as.vector(B)) %*% sdm %*% lambda +
                         diag(as.vector(B)) %*% apply(1 - sdm, 1, min), ncol = 1)

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

    return(list(raw.coeff = B, shrunk.coeff = matrix(B.shrunk, ncol = 1),
                lambda = lambda, Nrounds = nrounds, sdm = t(sdm)))
  }
}
