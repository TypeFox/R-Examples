#' Make Prediction
#'
#' Return risk prediction on new data.
#'
#' @param object [\code{rsig}]\cr
#'   An output object from rsig, see \code{\link[rsig]{rsig}}.
#' @param newdata [\code{data.frame}]\cr
#'   Data frame or matrix of input data (rows: examples, columns: features).
#' @param ... [ANY]\cr
#'   Additional arguments, currently ignored.
#' @return Risk prediction on new data.
#' @seealso \code{\link{rsig}}, \code{\link{rsig.eval}}, \code{\link{rsig.all}}
#' @export
#' @method predict rsig
#' @S3method predict rsig
predict.rsig = function(object, newdata, ...) {
  featnames = names(object$beta)
  newdata[, featnames, drop=FALSE] %*% object$beta + object$intercept
}

#' Performance Evaluation
#'
#' Evaluate performance on new data using predictions.
#'
#' @param pred [\code{predict.rsig}]\cr
#'   An output object from predict.rsig, see \code{\link[rsig]{predict.rsig}}.
#' @param surv.new [\code{Surv}]\cr
#'   Survival object, see \code{\link[survival]{Surv}}.
#' @param X.new [\code{data.frame}]\cr
#'   Data frame or matrix or matrix of input data (rows: examples, columns: features).
#' @param measures [\code{list}]\cr
#'   List of performance measures to be evaluated, "all" or in c("cindex", "tauc")
#' @param roc.time [\code{numeric(1)}]\cr
#'   Time to evaluate the time-dependent AUC. Default is \code{5}.
#' @seealso \code{\link{rsig}}, \code{\link{predict.rsig}}, \code{\link{rsig.all}}
#' @return Performance values
#' @export
rsig.eval = function(pred, surv.new, X.new, measures = "all", roc.time = 5) {
#   test = !train
  supported.measures = c("cindex", "tauc") #c("hr", "Dindex", "Cindex", "tROC", "brier")
  if (measures == "all")
    measures = supported.measures
  else
    measures = match.arg(measures, supported.measures, several.ok = TRUE)

  res = namedList(measures)
#   if ("hr" %in% measures) {
#     tmp = hazard.ratio(x = pred, surv.time = surv[test, 1L], surv.event = surv[test, 2L], na.rm = TRUE)
#     res$hr = list(hazard.ratio = tmp$hazard.ratio, pval = tmp$p.value)
#   }
#   if ("Dindex" %in% measures) {
#     tmp = D.index(x = pred, surv.time = surv[test, 1L], surv.event = surv[test, 2L], na.rm = TRUE)
#     res$Dindex = tmp$d.index
#   }
  if ("cindex" %in% measures) {
    tmp = concordance.index(x = pred, surv.time = surv.new[, 1L], surv.event = surv.new[, 2L], na.rm = TRUE, method = "noether")
    res$cindex = tmp$c.index
  }
  if ("tauc" %in% measures) {
    tmp = tdrocc(x = pred, surv.time = surv.new[, 1L], surv.event = surv.new[, 2L], na.rm = TRUE, time = roc.time)
#     res$troc = list(sens = tmp$sens, spec = tmp$spec, auc = tmp$AUC)
    res$tauc = tmp$AUC
  }
#   if ("brier" %in% measures) {
#     dd.tr = data.frame(time = surv[train, 1L], event = surv[test, 2L], score = pred)
#     dd.tt = data.frame(time = surv[test,  1L], event = surv[test, 2L], score = pred)
#     res$brier = sbrier.score2proba(data.tr = dd.tr, data.ts = dd.tt, method="cox")$bsc.integrated
#   }

  return(res)
}
