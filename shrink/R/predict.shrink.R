#' Predict Method for Objects of Class \code{shrink}
#'
#' Obtains predictions from shrunken regression coefficients from an object of
#' class \code{shrink}.
#' This class of objects is returned by the \code{shrink} function. Objects of this
#' class have methods for the functions \code{coef}, \code{predict}, \code{print},
#' \code{summary}, and \code{vcov}.

#' @param object an object of class \code{shrink}.
#' @param newdata a data frame for which predictions are obtained, otherwise
#'     predictions are based on the data stored in \code{object}.
#' @param type the type of prediction required.
#' @param shrinktype the type of shrinkage requested, if the \code{object} was
#'        obtained with \code{type = "all"}, either \code{"parameterwise"} or
#'        \code{"global"}.
#' @param terms with \code{type = "terms"} by default all terms are returned. A
#'     character vector specifies which terms are to be returned.
#' @param na.action function determining what should be done with missing values
#'     in \code{newdata}. The default is to include all observations.
#' @param collapse if \code{family = coxph} or \code{Cox}, an optional vector of
#'     subject identifiers. If specified, the output will contain one entry per
#'     subject rather than one entry per observation.
#' @param safe option from \code{predict.mfp}.
#' @param ... additional arguments to be passed to methods.
#'
#' @export
#' @importFrom utils getFromNamespace
#' @importFrom stats na.pass
#'
#' @note If \code{object} was obtained using \code{type = "all"}, \code{shrinktype}
#'       specifies for which type of shrinkage predictions are requested.
#'       \code{shrinktype} will be ignored if \code{object} was obtained using
#'       either \code{type = "parameterwise"} or \code{type = "global"}.
#'
#' @return A vector or matrix of predictions.
#'
#' @examples
#' data("GBSG")
#' library("mfp")
#'
#' fit <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05) +
#'            fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)
#'
#' dfbeta.global <- shrink(fit, type = "global",  method = "dfbeta")
#' dfbeta.pw     <- shrink(fit, type = "parameterwise", method = "dfbeta")
#' dfbeta.join   <- shrink(fit, type = "parameterwise", method = "dfbeta",
#'                         join=list(c("age.1", "age.2")))
#'
#' age <- 30:80
#' newdat <- data.frame(age = age, prm = 0)
#' refdat <- data.frame(age = 50, prm = 0)
#'
#' # unshrunken
#' plot(age, predict(fit, newdata = newdat, type = "lp") -
#'        predict(fit, newdata = refdat, type = "lp"), xlab = "Age",
#'      ylab = "Log hazard relative to 50 years", type = "l", lwd = 2)
#'
#' # globally shrunken
#' lines(age, predict(dfbeta.global,newdata = newdat, type = "lp") -
#'         predict(dfbeta.global, newdata = refdat, type = "lp"), lty = 3, col = "red", lwd = 2)
#'
#' # jointly shrunken
#' lines(age, predict(dfbeta.join, newdata = newdat, type = "lp") -
#'         predict(dfbeta.join, newdata = refdat, type = "lp"), lty = 4, col = "blue", lwd = 2)
#'
#' # parameterwise shrunken
#' lines(age, predict(dfbeta.pw, newdata = newdat, type = "lp") -
#'         predict(dfbeta.pw, newdata =refdat, type = "lp"), lty = 2, col = "green", lwd = 2)
#'
#' legend("topright", lty = c(1, 3, 4, 2), title = "SHRINKAGE",
#'        legend = c("No", "Global", "Joint", "Parameterwise"), inset = 0.01, bty = "n",
#'        col = c("black", "red", "blue", "green"), lwd = 2)
#' @seealso \code{\link{shrink}}, \code{\link{coef.shrink}}, \code{\link{print.shrink}},
#'     \code{\link{summary.shrink}}, \code{\link{vcov.shrink}}
predict.shrink <-
  function(object,
           newdata = NULL,
           type = c("link", "response", "lp", "risk", "expected", "terms"),
           shrinktype = NULL,
           terms = NULL,
           na.action = na.pass,
           collapse,
           safe = FALSE,
           ...)                             # se.fit = FALSE, dispersion = NULL
{
  # class glm = c("link", "response", "terms")
  # class mfp = c("link", "response", "lp", "risk", "expected", "terms"),
  # class coxph = c("lp", "risk", "expected", "terms")
  # class lm = c("response", "terms")

  if (!inherits(object, "shrink")) stop("'object' is not of class shrink")
  if (inherits(object$fit, "mfp")) { fit <- object$fit$fit } else { fit <- object$fit }

  if (object$type %in% "all") {
    if (is.null(shrinktype)) stop("Specify the type of shrinkage factors in 'shrinktype'; either 'global' or 'parameterwise'") else
    {
      shrinktype <- match.arg(shrinktype, c("global", "parameterwise"))
      if (shrinktype == "global") { fit$coefficients <- object$global$ShrunkenRegCoef } else
      { fit$coefficients <- object$parameterwise$ShrunkenRegCoef }
    }
  } else { fit$coefficients <- object$ShrunkenRegCoef }


  if (inherits(object$fit, "mfp")) {
    # code from mfp:predict.mfp Version 1.4.9 (Nov. 2012)
    if (object$fit$family$family == "Cox") {
      if (is.null(object$fit$terms)) terms = names(object$fit$assign)
      if (!missing(newdata))
        if (!missing(collapse))
          getFromNamespace("predict.coxph", "survival")(fit, newdata = newdata, type = type,
                                                        se.fit = FALSE, terms = terms,
                                                        collapse = collapse, safe = safe, ...)
      else
        getFromNamespace("predict.coxph", "survival")(fit, newdata = newdata, type = type,
                                                      se.fit = FALSE, terms = terms,
                                                      safe = safe, ...)
      else
        if (!missing(collapse))
          getFromNamespace("predict.coxph", "survival")(fit, type = type, se.fit = FALSE,
                                                        terms = terms, collapse = collapse,
                                                        safe = safe, ...)
      else
        getFromNamespace("predict.coxph", "survival")(fit, type = type, se.fit = FALSE,
                                                      terms = terms, safe = safe, ...)
    } else {
      if (!missing(newdata))
        predict.glm(fit, newdata = newdata, type = type, se.fit = FALSE,
                    dispersion = NULL, terms = terms, na.action = na.action, ...)
      else
        predict.glm(fit, type = type, se.fit = FALSE, dispersion = NULL,
                    terms = terms, na.action = na.action, ...)
    }
  } else {
    if (inherits(object$fit, c("glm", "lm"))) {
      if (!missing(newdata)) {
        predict.glm(fit, newdata = newdata, type = type, se.fit = FALSE,
                    dispersion = NULL, terms = terms, na.action = na.action, ...)
      } else {
        predict.glm(fit, type = type, se.fit = FALSE, dispersion = NULL,
                    terms = terms, na.action = na.action, ...)
      }
    } else {
      #    if (inherits(object$fit, "coxph")) {
      if (!is.null(terms)) { terms <- names(fit$assign) }
      if (is.null(newdata)) { newdata <- data.frame(fit$x) }
      predict(fit, newdata = newdata, type = type, se.fit = FALSE,
              na.action = na.action, terms = terms, collapse,  ...)
    }
  }
}
