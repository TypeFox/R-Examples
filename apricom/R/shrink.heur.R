#' Shrinkage After Estimation Using Heuristic Formulae
#'
#' Shrink regression coefficients using heuristic formulae, first described
#' by Van Houwelingen and Le Cessie (Stat. Med., 1990)
#'
#' This function can be used to estimate shrunken regression coefficients based on
#' heuristic formulae (see References). A linear or logistic regression model with
#' an intercept is fitted to the data, and a shrinkage factor is estimated. The
#' shrinkage factor is then applied to the regression coefficients. If
#' \code{int.adj == FALSE} the intercept value is estimated as described in
#' Harrell 2001.If \code{int.adj == TRUE} the intercept value will be re-estimated
#' by refitting the model with the shrunken coefficients.
#'
#' The heuristic formula work by applying the
#' number of model degrees of freedom (or the number of predictors) as a penalty,
#' and so as the model becomes more complex, the necessary shrinkage increases, and
#' the shrinkage factor becomes closer to zero.
#'
#' @importFrom stats glm.fit binomial
#'
#' @param dataset   a dataset for regression analysis. Data should be in the form
#' of a matrix, with the outcome variable as the final column. Application of the
#' \code{\link{datashape}} function beforehand is recommended, especially if
#'  categorical predictors are present. For regression with an intercept included
#'  a column vector of 1s should be included before the dataset (see examples)
#' @param model     type of regression model. Either "linear" or "logistic".
#' @param DF        the number of degrees of freedom or number of predictors in the model.
#'                  If DF is missing the value will be automatically estimated. This
#'                  may be inaccurate for complex models with non-linear terms.
#' @param int       logical. If TRUE the model will include a regression intercept.
#' @param int.adj   logical. If TRUE the regression intercept will be re-estimated
#'          after shrinkage of the regression coefficients. If FALSE the regression
#'          intercept will be re-estimated as described by Harrell 2001.
#'
#' @return \code{shrink.heur} returns a list containing the following:
#' @return \item{raw.coeff}{the raw regression model coefficients, pre-shrinkage.}
#' @return \item{shrunk.coeff}{the shrunken regression model coefficients}
#' @return \item{lambda}{the heuristic shrinkage factor}
#' @return \item{DF}{the number of degrees of freedom or number of predictors in the model}
#'
#' @note Warning: In poorly fitting models that includea large number of predictors
#'       the number of degrees of freedom may approch or exceed the model chi square.
#'       In such cases the shrinkage factor will be very small or even negative,
#'       and a different model building strategy is recommended.
#'
#' @examples
#'## Example 1: Linear regression using the iris dataset
#'## shrinkage using a heuristic formula
#' data(iris)
#' iris.data <- as.matrix(iris[, 1:4])
#' iris.data <- cbind(1, iris.data)
#' set.seed(123)
#' shrink.heur(dataset = iris.data, model = "linear")
#'
#'## Example 2: logistic regression using a subset of the mtcars data
#'## shrinkage using a heuristic formula
#' data(mtcars)
#' mtc.data <- cbind(1,datashape(mtcars, y = 8, x = c(1,6,9)))
#' head(mtc.data)
#' set.seed(321)
#' shrink.heur(dataset = mtc.data, model = "logistic", DF = 3,
#' int = TRUE, int.adj = TRUE)
#'
#' @references Harrell, F. E. \emph{"Regression modeling strategies: with applications
#'              to linear models, logistic regression, and survival analysis."} \emph{Springer}, (2001).
#' @references Harrell, F. E., Kerry L. Lee, and Daniel B. Mark. \emph{"Tutorial in
#'            biostatistics multivariable prognostic models: issues in developing models,
#'            evaluating assumptions and adequacy, and measuring and reducing errors."}
#'            \emph{Statistics in medicine} (1996) 15:361-387.
#' @references Steyerberg, E. \emph{"Clinical Prediction Models"} \emph{Springer} (2009)
#' @references Van Houwelingen, J. C. and Le Cessie, S., \emph{"Predictive value of statistical models."}
#'             \emph{Statistics in medicine} (1990) 9:1303:1325.

shrink.heur <- function(dataset, model, DF, int = TRUE, int.adj = FALSE){

  if(missing(model)) stop("Type of model must be specified in the function call")
  if(int == FALSE) stop("For heuristic methods, an intercept must be included in
                        the regression")
  if(missing(DF)) warning("The number of degrees of freedom used by predictors was
                          estimated automatically")
  if (missing(DF)) DF <- dim(dataset)[2] - 2

  dataset <- as.matrix(dataset)
  nc <- dim(dataset)[2]


  if (model == "linear") {

    ### HEURISTIC FORMULA FOR LINEAR REGRESSION WITH OLS ###

    Rsq <- function(dataset){
      nc <- dim(dataset)[2]
      b <- ols.rgr(dataset)
      y <- dataset[, nc]
      X <- dataset[, 1:(nc - 1)]
      yhat <- (X %*% b)
      res <- y - yhat
      sse0 <- t(res) %*% res
      scale(y, center = TRUE, scale = FALSE)
      Syy <- t(y) %*% y
      r <- 1 - (sse0 / Syy)
    }
    b <- ols.rgr(dataset)
    R2<- Rsq(dataset)
    n <- dim(dataset)[1]
    AdjR2 <- 1 - (1 - R2) * (n - 1) / (n - DF - 1)
    s <- (n - DF - 1) / (n - 1) * AdjR2 / R2

  } else if (model == "logistic") {

    ### HEURISTIC FORMULA FOR GLMS WITH MAXIMUM LIKELIHOOD ESTIMATION ###

    # calculate the model coefficients using MLE
    b <- ml.rgr(dataset)
    # do the same but with an empty model (intercept only)
    y <- dataset[, nc]
    d0 <- cbind(1, y)
    b0 <- ml.rgr(d0)
    b0 <- matrix(c(b0, rep(0, nc - 2)), ncol = 1)
    model.chisq <- loglikelihood(b0, dataset) - loglikelihood(b, dataset)
    s <- ((model.chisq - DF) / model.chisq)
  }

  ### APPLY SHRINKAGE FACTOR ### As described in Harrell 2001 (book p 64)
  if (int.adj == FALSE) {
    b.shrunk <- b
    y <- dataset[, nc]
    bl <- dim(b.shrunk)[1]
    b.shrunk[1] <- (1 - s) * mean(y) + s * b[1]
    b.shrunk[2:bl] <- b[2:bl] * s

    } else {
      b.shrunk <- matrix(s * b[2:length(b)], ncol = 1)
      if (model == "linear") {
        new.int <- mean((dataset[, nc]) - ((dataset[, 2:(nc - 1)]) %*% b.shrunk))
        b.shrunk <- c(new.int, b.shrunk)
      } else {
        X.i <-  matrix(dataset[, 1], ncol = 1)
        Y <- dataset[, nc]
        offs <- as.vector((dataset[, 2:(nc - 1)]) %*% b.shrunk)
        new.int <- glm.fit(X.i, Y, family = binomial(link = "logit"),
                           offset = offs)$coefficients
        b.shrunk <- c(new.int, b.shrunk)
        }
    }

  return(list(raw.coeff = b, shrunk.coeff = matrix(b.shrunk, ncol = 1),
              lambda = s, DF = DF))
}

