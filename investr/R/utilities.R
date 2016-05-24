#' 'Not Available' / Missing Values
#' 
#' A simple version of the newer \code{anyNA} function which essentially
#' implements any(is.na(x)).
#' 
#' @keywords internal
AnyNA <- function(x) {
  any(is.na(x))
}
  

#' Make new data frame
#' 
#' Create a new data frame from a specified x value that has the same structure 
#' as the data frame used to create \code{object}. (For internal use only.)
#' @importFrom stats setNames
#' @keywords internal
makeData <- function(x, label) {
  setNames(data.frame(x), label)
}


#' Design Matrix for Fixed and Random Effects
#'
#' Construct the fixed or random effects design matrix from a fitted model using 
#' \code{newdata}. (For internal use only.)
#' 
#' @keywords internal
makeX <- function(object, newdata) {
  model.matrix(eval(object$call$fixed)[-2], data = newdata)
}


#' @rdname makeX
#' @importFrom nlme asOneFormula
#' @importFrom stats formula na.fail model.matrix
#' @keywords internal
makeZ <- function(object, newdata) {
  Q <- object$dims$Q  # number of grouping levels
  mCall <- object$call  # list containing image of the nlme call
  fixed <- eval(eval(mCall$fixed)[-2])  # fixed effects formula
  reSt <- object$modelStruct$reStruct  # random effects structure
  mfArgs <- list(formula = asOneFormula(formula(reSt), fixed),
                 data = newdata, na.action = na.fail,
                 drop.unused.levels = TRUE)
  dataMix <- do.call("model.frame", mfArgs)
  model.matrix(reSt, dataMix)
}


#' Extract residual standard error
#' 
#' Extract residual standard error from a fitted model. (For internal use 
#' only.)
#' 
#' @keywords internal
Sigma <- function(object, ...) {
  UseMethod("Sigma")
} 
Sigma.lm <- function(object, ...) summary(object)$sigma
Sigma.nls <- function(object, ...) summary(object)$sigma
Sigma.lme <- function(object, ...) object$sigma


#' Evaluate response variance
#'
#' Evaluate response variance at a given value of the predictor variable. (For 
#' internal use only.)
#' 
#' @importFrom nlme getVarCov
#' @keywords internal
varY <- function(object, newdata) {
  Zmat <- makeZ(object, newdata)  # random effects design matrix
  Gmat <- getVarCov(object)  # random effects variance-covariance matrix
  var.y <- Zmat %*% Gmat %*% t(Zmat) + Sigma(object)^2  # ZGZ' + (sigma^2)I
  if (is.matrix(var.y)) unname(diag(var.y)) else var.y
}