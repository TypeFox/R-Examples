#' Predicted Values from dosresmeta Models 
#'
#' @description This method function computes predictions from fitted dose-response models 
#' represented in objects of class "\code{dosresmeta}", optionally for a new set of exposure levels. 
#' Predictions are optionally accompanied by confidence intervals and/or standard errors for the predictions.
#' 
#' @param object an object of class \code{dosreseta}.
#' @param newdata an optional data frame or matrix in which to look for variables values with which to predict from dose-response models.
#' @param xref an optional scalar to indicate which levels should serve as referent for the predicted relative risks. See details.
#' @param expo logical switch indicating if the prediction should be on the exponential scale.
#' @param se.incl logical switch indicating if standard errors need to be included.
#' @param ci.incl logical switch indicating if confidence intervals need to be included.
#' @param ci.level a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.
#' @param order logical to indicate if the predictions need to be sorted by exposure levels.
#' @param delta an optional scalar to specify to predict the linear trend related to that increase.
#' @param \dots further arguments passed to or from other methods.
#' @details The method function \code{predict} produces predicted values from \code{dosresmeta} objects. When more than one study is included in the analysis,
#' estimated predictions are only based on the fixed part of the model.
#' 
#' If \code{newdata} is omitted, the predictions are based on the data used for the fit. If \code{xref} is provided, it must be equal to one of the modeled values. If not
#' provided, the minimum modeled referent value will be used as referent for the predicted relative risks
#' 
#' If \code{newdata} is specified, it should include all the variables used to model the dose-response relation. Again, if specified, \code{xref} must be equal to one
#' of the value in the newdata. If omitted, the minimum value for the newdara will be used as referent.
#' 
#' Only for the linear trend it is possible to specify the predicted increase of risk correspongind to an increase equal to \code{delta} argument.
#' 
#' By default (\code{order = TRUE}), the predictions are sorted by exposure levels to facilitate understanding and possible graphical
#' presentation of the results.
#' 
#' @return The results are returned structured in a data frame.
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' @seealso \code{\link{dosresmeta}}, \code{\link{predict}}
#' 
#' @examples
#' ## Load data and run the model
#'data("alcohol_cvd")
#'model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                    se = se, cases = cases, n = n, data = alcohol_cvd) 
#'                    
#'## Predicted modeled data
#'predict(model, order = FALSE)
#'
#'## Plot predicted dose-response relation
#'with(predict(model), {
#'  plot(dose, pred, log = "y", type = "l",
#'       xlim = c(0, 45), ylim = c(.4, 2))
#'  lines(dose,  ci.lb, lty = 2)
#'  lines(dose, ci.ub, lty = 2)
#'  rug(dose, quiet = TRUE)
#'})
#'
#'## Prediction for new values
#'newdata <- data.frame(dose = seq(0, 50, 1))
#'predict(model, newdata)
#'
#'## Smoother plot
#'with(predict(model, newdata),{
#'  plot(dose, pred, log = "y", type = "l",
#'                ylim = c(.4, 2))
#'  lines(dose, ci.lb, lty = 2)
#'  lines(dose, ci.ub, lty = 2)
#'  rug(alcohol_cvd$dose, quiet = TRUE)
#'})
#'
#'## Tabular results
#'newdata <- data.frame(dose = seq(0,50,5))
#'round(predict(model, newdata), 2)
#' 
#' @rdname predict.dosresmeta
#' @method predict dosresmeta
#' @S3method predict dosresmeta
#' @export predict.dosresmeta
#' 
predict.dosresmeta <- function(object, newdata, xref, se.incl = FALSE, expo = TRUE, 
                               ci.incl = TRUE, ci.level = .95, order = TRUE, delta, ...)
{

  if (!missing(delta)){
    newdata <- data.frame(c(0, delta))
    colnames(newdata) <- colnames(object$coefficients)
  }
  if (missing(newdata) || is.null(newdata)) {
    X <- as.matrix(object$design)
    if (missing(xref)){
      xref <- data.frame(object$xref)
    }
  } else{
    Terms = delete.response(object$termsi)
    m <- model.frame(Terms, newdata)
    X <- model.matrix(Terms, m)
    if (missing(xref)){
      xref <- min(X[, 1])
    }
  }
  X0 <- matrix(X[X[, 1] == as.double(xref), ], ncol = ncol(X))[1,]    
  Xref <- X - matrix(rep(X0, nrow(X)), ncol = ncol(X), byrow = T) 
  pred <- tcrossprod(Xref, object$coefficients)
  fit <- if (expo == T){
    cbind(X, exp(pred))
  } else cbind(X, pred)
  colnames(fit) <- c(colnames(X), "pred")
  zvalci <- qnorm((1 - ci.level)/2, lower.tail = FALSE)
  if (ci.incl == T){
    se <- diag(Xref %*% vcov(object) %*% t(Xref))^.5
    ci.lb <- pred -  zvalci * se
    ci.ub <- pred +  zvalci * se
    fit <- if (expo == T){
      cbind(fit, exp(ci.lb), exp(ci.ub))
    } else cbind(fit, ci.lb, ci.ub)
    colnames(fit) <- c(colnames(X), "pred", "ci.lb", "ci.ub")
  }
  if (se.incl == T) fit <- cbind(fit, se = se)
  fit <- as.data.frame(fit)
  if (order == T) fit <- fit[order(fit[, 1]),]
  if (!missing(delta)){
    fit <- fit[2, ]
    rownames(fit) <- ""
  }
  fit
}