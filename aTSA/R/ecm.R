#' Error Correction Model
#' @description Fits an error correction model for univriate response.
#' @param y a response of a numeric vector or univariate time series.
#' @param X an exogenous input of a numeric vector or a matrix for multivariate time series.
#' @param output a logical value indicating to print the results in R console. 
#' The default is \code{TRUE}.
#' @details An error correction model captures the short term relationship between the 
#' response \code{y} and the exogenous input variable \code{X}. The model is defined as 
#' \deqn{dy[t] = bold{\beta}[0]*dX[t] + \beta[1]*ECM[t-1] + e[t],}
#' where \eqn{d} is an operator of the first order difference, i.e., 
#' \eqn{dy[t] = y[t] - y[t-1]}, and \eqn{bold{\beta}[0]} is a coefficient vector with the 
#' number of elements being the number of columns of \code{X} (i.e., the number
#' of exogenous input variables), and\eqn{ ECM[t-1] = y[t-1] - hat{y}[t-1]} which is the 
#' main term in the sense that its coefficient \eqn{\beta[1]} explains the short term 
#' dynamic relationship between \code{y} and \code{X} 
#' in this model, in which \eqn{hat{y}[t]} is estimated from the linear regression model 
#' \eqn{y[t] = bold{\alpha}*X[t] + u[t]}. Here, \eqn{e[t]} and \eqn{u[t]} are both error terms
#' but from different linear models.
#' 
#' @note Missing values are removed before the analysis. In the results, \code{dX} or 
#' \code{dX1}, \code{dX2}, ... represents the first difference of each exogenous input 
#' variable \code{X}, and \code{dy} is the first difference of response \code{y}.
#' @return An object with class "\code{lm}", which is the same results of \code{\link{lm}} for 
#' fitting linear regression.
#' @author Debin Qiu
#' @references
#' Engle, Robert F.; Granger, Clive W. J. (1987). Co-integration and error correction: 
#' Representation, estimation and testing. \emph{Econometrica}, 55 (2): 251-276.
#' @examples X <- matrix(rnorm(200),100,2)
#' y <- 0.1*X[,1] + 2*X[,2] + rnorm(100)
#' ecm(y,X)
#' @importFrom stats lm
#' @importFrom stats residuals
#' @export

ecm <- function(y,X,output = TRUE)
{
  if (NCOL(y) > 1)
    stop("'y' must be a numeric vector or univariate time series")
  if (any(!is.finite(y)) || any(!is.finite(X)))
    stop("missing values are not allowed for 'y' and 'X'")
  n.y <- length(y)
  n.X <- NROW(X)
  m.X <- NCOL(X)
  if (n.y != n.X)
    stop("'y' and 'X' must have same number of observations")
  if (n.y < 2L)
    stop("'y' must have at least two observations")
  dy <- diff(y,1)
  dX <- if (m.X > 1) apply(X,2,"diff",1) else diff(X,1)
  res <- residuals(lm(y ~ X - 1))
  ECM <- res[-1]
  lm.fit <- lm(dy ~ dX + ECM - 1)
  lm.sum <- summary(lm.fit)
  lm.sum$coefficients[m.X + 1,1] <- -lm.sum$coefficients[m.X + 1,1]
  coefs <- lm.sum$coefficients
  if (output) print(lm.sum)
  ecm <- lm.sum
}