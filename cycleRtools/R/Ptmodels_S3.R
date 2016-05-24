# S3 Methods
# ------------------------------------------------------------------------------
#' @export
print.Ptmodels <- function(x, ...) {
  print.data.frame(x$table, digits = 2)
}
#' @export
coef.Ptmodels <- function(object, ...) {
  lapply(object$models, coef)
}
#' @export
summary.Ptmodels <- function(object, ...) {
  lapply(object$models, summary)
}
#' Predict Power or Time
#'
#' Given a Ptmodels \code{object}, the predict.Ptmodels will produce a named
#' numeric vector of either time (seconds) or power (watts) values according to
#' the \code{x} and \code{y} arguments
#'
#' @param object an object of class "Ptmodels".
#' @param x the value for which to make a prediction.
#' @param xtype what is \code{x}? A power or a time value?
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#' data(Pt_prof)  # Example power-time profile.
#'
#' P    <- Pt_prof$pwr
#' tsec <- Pt_prof$time
#'
#' mdls <- Pt_model(P, tsec)  ## Model.
#' print(mdls)
#'
#' ## What is the best predicted 20 minute power?
#' predict(mdls, x = 60 * 20, xtype = "time")
#'
#' ## How sustainable is 500 Watts?
#' predict(mdls, x = 500, xtype = "P") / 60  # Minutes.
#'
#' ## Create some plots of the models.
#' par(mfrow = c(2, 2), mar = c(3.1, 3.1, 1.1, 1.1))
#' plotargs <- alist(x = tsec, y = P, cex = 0.2, ann = FALSE, bty = "l")
#' mapply(function(f, m) {
#'   do.call(plot, plotargs)
#'   curve(f(x), col = "red", add = TRUE)
#'   title(main = paste0(rownames(m),"; RSE = ", round(m$RSE, 2)))
#'   legend("topleft", legend = m$formula, bty = "n")
#'   return()
#' }, f = mdls$Pfn, m = split(mdls$table, seq_len(nrow(mdls$table))))
#'
#' @return a named numeric vector of predicted values. Names correspond to their
#'   respective models.
#'
#' @export
predict.Ptmodels <- function(object, x, xtype = c("pwr", "time"), ...) {
  if (length(xtype) > 1) xytpe <- xtype[1]
  fn <- switch(xtype,  # Match multiple inputs.
               "P" = "tfn", "power" = "tfn", "pwr" = "tfn",
               "tsec" = "Pfn", "time" = "Pfn", "t" = "Pfn",
               "Pfn")  # Default.
  out <- suppressWarnings(vapply(object[[fn]], function(f) f(x), numeric(1)))
  if (fn == "tfn") out[out < 0 | is.nan(out)] <- Inf
  out
}
