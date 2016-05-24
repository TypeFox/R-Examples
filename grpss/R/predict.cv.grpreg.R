#' Model predictions for a fitted \code{cv.grpreg} object
#' @description Similar to usual predict methods and
#' \code{predict.cv.grpreg} in \code{grpreg} package.
#' @param object A fitted "\code{grpreg}" object from \code{\link{grpss}},
#' or \code{\link[grpreg]{cv.grpreg}} function.
#' @param newdata Optionally, a matrix or data frame where to predict. If omits, the fitted
#' predictors are used.
#' @param type The type of prediction: "\code{response}" gives the fitted values; \code{"class"}
#' returns the predicted class for the binomial outcome; "\code{probability}" returns the
#' predicted probabilities for the logistic regression.
#' @param ... Not used.
#' @details This function gives the predictions at \code{newdata} or all predictors if the
#' argument \code{newdata} is not supplied. Typically, \code{type = "response"} is
#' used for linear or poisson regression, and \code{type = "class"} or
#' \code{type = "probability"} is used for logistic regression.
#' @return The predicted values depending on the type.
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grpss}}
#' @examples
#' library(MASS)
#' set.seed(23)
#' n <- 30 # sample size
#' p <- 3  # number of predictors in each group
#' J <- 50  # group size
#' group <- rep(1:J,each = 3)  # group indices
#' X <- mvrnorm(n,seq(0,5,length.out = p*J),diag(p*J))
#' beta <- runif(12,-2,5)
#' mu <- X%*%matrix(c(beta,rep(0,p*J-12)),ncol = 1)
#'
#' # linear regression with family = "gaussian"
#' y <-  mu + rnorm(n)
#'
#' ## with cross-validation
#' gss11 <- grpss(X,y,group,select = TRUE,cross.validation = TRUE)
#' predict(gss11)  # fitted values
#' predict(gss11, newdata = X[1,])  # predicted values at given 'newdata'
#'
#' # logistic regression with family = "binomial"
#' set.seed(23)
#' y1 <- rbinom(n,1,1/(1 + exp(-mu)))
#' gss21 <- grpss(X,y1,group, criterion = "gDC",select = TRUE,
#'                family = "binomial")
#' predict(gss21)
#'
#' @export

predict.cv.grpreg <- function(object,newdata,
                              type = c("response","class","probability"), ...) {
  predict.grpreg(object,newdata,lambda = NULL,type,...)
}
