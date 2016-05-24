#' Calculate performance for paired bootstrapped ROC curves
#' 
#' For a given metric this calculates the difference in performance between two paired predictors
#' stored in an object of class \code{fbroc.paired.roc} in addition to their individual performance.
#' 
#' @param roc An object of class \code{fbroc.paired.roc}.
#' @inheritParams perf.fbroc.roc
#' @param metric A performance metric. Select "auc" for the AUC, "tpr" for the TPR at a fixed
#' FPR and "fpr" for the FPR at a fixed TPR.
#' @export
#' @examples
#' data(roc.examples)
#' example <- boot.paired.roc(roc.examples$Cont.Pred, roc.examples$Cont.Pred.Outlier,
#'                                roc.examples$True.Class)
#' perf(example, metric = "auc")   
#' perf(example, metric = "tpr", fpr = 0.2) # Get difference in TPR at a FPR of 20%                             
perf.fbroc.paired.roc <- function(roc, metric = "auc", conf.level = 0.95, tpr = NULL, fpr = NULL, ...) {
#perf.paired.roc <- function(roc, metric = "auc", conf.level = 0.95, tpr = NULL, fpr = NULL) {
  # start with data validation
  if (!is(roc, "fbroc.paired.roc"))
    stop("roc must be of class fbroc.paired.roc")
  if (length(metric) != 1 | class(metric) != "character")
    stop("metric must be character")
  if (!(metric %in% c("auc", "tpr", "fpr")))
    stop(paste(metric,"is not a valid performance metric"))
  
  if (metric == "auc") {
    metric.text = "AUC"
    metric.number <- as.integer(0)
    param.vec <- 0
  }
  if (metric == "tpr") {
    param.vec = validate.single.numeric(fpr, "FPR")
    metric.text <- paste("TPR at a fixed FPR of", round(param.vec, 3))
    metric.number <- as.integer(1)
  }
  if (metric == "fpr") {
    param.vec = validate.single.numeric(tpr, "TPR")
    metric.text <- paste("FPR at a fixed TPR of", round(param.vec, 3))
    metric.number <- as.integer(2)
  }
  
  
  # call C++ to calculate actual results
  tpr.m1 <- matrix(roc$roc1$TPR, nrow = 1)
  fpr.m1 <- matrix(roc$roc1$FPR, nrow = 1)
  tpr.m2 <- matrix(roc$roc2$TPR, nrow = 1)
  fpr.m2 <- matrix(roc$roc2$FPR, nrow = 1)
  
  observed.perf1 <- get_cached_perf(tpr.m1, fpr.m1, param.vec, metric.number)
  observed.perf2 <- get_cached_perf(tpr.m2, fpr.m2, param.vec, metric.number)
  if (roc$use.cache) {
    perf.boot.list <- get_cached_perf_paired(roc$boot.tpr1, roc$boot.fpr1, 
                                             roc$boot.tpr2, roc$boot.fpr2,
                                             param.vec, metric.number)
  } else {
    perf.boot.list <- get_uncached_perf_paired(roc$predictions1,
                                               roc$predictions2,
                                               as.integer(roc$true.classes),
                                               param.vec,
                                               roc$n.boot,
                                               metric.number)
  }
  
  # Quantile based confidence interval
  alpha <- 0.5 * (1 - conf.level)
  alpha.levels <- c(alpha, 1 - alpha) 
  ci1 <- as.numeric(quantile(perf.boot.list[[1]], alpha.levels))
  names(ci1) <- NULL
  ci2 <- as.numeric(quantile(perf.boot.list[[2]], alpha.levels))
  names(ci2) <- NULL
  ci.diff <- as.numeric(quantile(perf.boot.list[[1]] - perf.boot.list[[2]], alpha.levels))
  names(ci.diff) <- NULL
  cor.pred <- cor(perf.boot.list[[1]], perf.boot.list[[2]])
  
  perf <- list(Observed.Performance.Predictor1 = observed.perf1,
               CI.Performance.Predictor1 = ci1,
               Observed.Performance.Predictor2 = observed.perf2,
               CI.Performance.Predictor2 = ci2,
               Observed.Difference = observed.perf1 - observed.perf2,
               CI.Performance.Difference = ci.diff,
               conf.level = conf.level,
               Cor = cor.pred,
               metric = metric.text,
               params = param.vec,
               n.boot = roc$n.boot,
               boot.results.pred1 = perf.boot.list[[1]],
               boot.results.pred2 = perf.boot.list[[2]]
  )
  class(perf) <- append(class(perf), "fbroc.perf.paired")
  return(perf)
}