#' Bootstrap ROC curve
#'
#' \code{boot.roc} calculates the ROC curve, initializes the settings
#' and calculates the bootstrap results for the true and false
#' positive rate at every relevant threshold. Missing values are removed with 
#' a warning prior to bootstrapping.
#'
#' @param pred A numeric vector. Contains predictions. \code{boot.roc} assumes
#'   that a high prediction is evidence for the observation belonging to the
#'   positive class.
#' @param true.class A logical vector. TRUE indicates the sample belonging to the
#'   positive class.
#' @param stratify Logical. Indicates whether we use stratified bootstrap.
#'   Default to TRUE. Non-stratified bootstrap is not yet implemented.
#' @param n.boot A number that will be coerced to integer. Specified the 
#'   number of bootstrap replicates. Defaults to 1000.
#' @param use.cache If true the bootstrapping results for the
#'   ROC curve will be pre-cached. This increases speed when the object is used often, but also
#'   takes up more memory.
#' @param tie.strategy How to handle ties. See details below.
#' @return A list of class \code{fbroc.roc}, containing the elements:
#' \item{prediction}{Input predictions.}
#' \item{true.class}{Input classes.}
#' \item{roc}{A data.frame containing the thresholds of the ROC curve and the TPR and FPR at these
#' thresholds.}
#' \item{n.thresholds}{Number of thresholds.}
#' \item{n.boot}{Number of bootstrap replicates.}
#' \item{use.cache}{Indicates if cache is used for this ROC object}
#' \item{tie.strategy}{Used setting how to handle ties in predictors.}
#' \item{n.pos}{Number of positive observations.}
#' \item{n.neg}{Number of negative observations.}
#' \item{auc}{The AUC of the original ROC curve.}
#' \item{boot.tpr}{If the cache is enabled, a matrix containing the bootstrapped TPR at the thresholds.}
#' \item{boot.fpr}{If the cache is enabled, a matrix containing the bootstrapped FPR at the thresholds.}
#' @section Caching:
#' If you enable caching, \code{boot.roc} calculates the requested number of bootstrap samples and
#' saves the TPR and FPR values for each iteration. This can take up a sizable portion of memory,
#' but it speeds up subsequent operations. This can be useful if you plan to use the ROC curve
#' multiple \code{fbroc} functions.
#' @section Ties:
#' You can set this parameter to either 1 or 2. If your numerical predictor has no ties, both settings
#' will produce the same results. 
#' If you set \code{tie.strategy} to 1 the ROC curve is built by connecting the TPR/FPR pairs for
#' neighboring thresholds. A tie.strategy of 2 indicates that the TPR calculated at a specific FPR
#' is the best TPR at a FPR smaller than or equal than the FPR specified. Defaults to 2.
#' @examples
#' y <- rep(c(TRUE, FALSE), each = 500)
#' x <- rnorm(1000) + y
#' result.boot <- boot.roc(x, y)
#' @seealso \url{http://www.epeter-stats.de/roc-curves-and-ties/}, \code{\link{plot.fbroc.roc}}, 
#' \code{\link{print.fbroc.roc}}
#' 
#' @export
boot.roc <- function(pred, true.class, stratify = TRUE, n.boot = 1000,
                     use.cache = FALSE, tie.strategy = NULL) {
  # validate input
  if ((length(pred) != length(true.class)))
    stop("Predictions and true classes need to have the same length")
  if (class(pred) == "integer") {
    pred <- as.numeric(pred)
  }
  if ((class(pred) != "numeric"))
    stop("Predictions must be numeric")
  if ((class(true.class) != "logical"))
    stop("Classes must be logical")
  if ((class(stratify) != "logical"))
    stop("Classes must be logical")
  
  index.na <- is.na(pred) | is.na(true.class)
  if (any(index.na)) {
    n <- sum(index.na)
    warning.msg <- 
      paste(n, "observations had to be removed due to missing values")
    warning(warning.msg)
    true.class <- true.class[!index.na]
    pred <- pred[!index.na]
  }
  if (is.null(tie.strategy)) {
    if (length(pred) == length(unique(pred))) tie.strategy <- 1 
    else tie.strategy <- 2
  }
  
  if (!(tie.strategy %in% 1:2)) stop("tie.strategy must be 1 or 2")
  if (sum(true.class) == 0)
    stop("No positive observations are included")
  if (sum(!true.class) == 0)
    stop("No negative observations are included")
  
  n.boot <- as.integer(n.boot)
  
  if (length(n.boot) != 1)
    stop("n.boot must have length 1")
  
  if (length(stratify) != 1)
    stop("stratify must have length 1")
  
  if (!stratify) stop("Non-stratified bootstrapping is not yet supported")
  
  true.int <- as.integer(true.class)

  original.roc <- roc_analysis(pred, true.int)
  auc <- original.roc[[4]]
  original.roc <- c.list.to.roc(original.roc)
  if (use.cache) {
    booted.roc <- tpr_fpr_boot(pred, true.int, n.boot)
    boot.tpr <- booted.roc[[1]]
    boot.fpr <- booted.roc[[2]]
    rm(booted.roc)
  } else {
    boot.tpr <- NULL
    boot.fpr <- NULL
  }
  
  output <- list(predictions = pred,
                 true.classes = true.class,
                 n.thresholds = nrow(original.roc),
                 n.boot = as.integer(n.boot),
                 use.cache = use.cache,
                 tie.strategy = tie.strategy,
                 n.pos = sum(true.class),
                 n.neg = sum(!true.class),
                 roc = original.roc,
                 auc = auc,
                 boot.tpr = boot.tpr,
                 boot.fpr = boot.fpr)
  class(output) <- append(class(output), "fbroc.roc")
  return(output)
}


conf.roc <- function(roc, conf.level = 0.95, steps = 250) {
  .Deprecated("conf")
  perf(conf, conf.level = conf.level, steps = steps)
}

#' Generates confidence intervals for the TPR for a range of FPRs or vice versa
#' 
#' Calculates confidence intervals for the TPR at different FPR values or vice versa. The stepsize
#' at which the TPR or FPR is calculated can be set as needed.
#' 
#' @param roc Object of class \code{fbroc.roc}.
#' @param conf.level Confidence level to be used for the confidence intervals. Defaults to 0.95.
#' @param conf.for Use "tpr" to get confidence regions for the TPR at specific FPRs. Use "fpr"
#' instead for confidence regions for the FPR at specific TPRs.
#' @param steps Number of discrete steps at which the requested rate and the confidence region is calculated.
#' Defaults to 250.
#' @param ... Further arguments, that are not used at this time.
#' @return A data.frame containing either discrete TPR steps and estimates and confidence bounds for
#' FPR or vice versa, depending upon \code{conf.for}.
#' @export
#' @examples 
#' data(roc.examples)
#' example <- boot.roc(roc.examples$Cont.Pred, roc.examples$True.Class,
#'                     n.boot = 100)
#' conf(example, conf.for = "tpr", steps = 10) # get confidence regions for TPR at FPR
#' conf(example, conf.for = "fpr", steps = 10) # get confidence regions for FPR at TPR
#' @seealso \code{\link{boot.roc}}
conf.fbroc.roc <- function(roc, conf.level = 0.95, conf.for = "tpr", steps = 250, ...) {
  if (!(conf.for %in% c("tpr", "fpr"))) stop("Invalid rate given for confidence region")
  alpha <- 0.5*(1 - conf.level)
  alpha.levels <- c(alpha, 1 - alpha) 
  steps = as.integer(steps)
  
  if (conf.for == "tpr") {
    cached.fun <- tpr_at_fpr_cached
    uncached.fun <- tpr_at_fpr_uncached
    frame.names <- c("FPR", "TPR", "Lower.TPR", "Upper.TPR")
  } else {
    cached.fun <- fpr_at_tpr_cached
    uncached.fun <- fpr_at_tpr_uncached
    frame.names <- c("TPR", "FPR", "Lower.FPR", "Upper.FPR")
  }
  estimate <- cached.fun(matrix(roc$roc$TPR, nrow = 1), 
                         matrix(roc$roc$FPR, nrow = 1), 
                         steps)
  if (roc$use.cache) {
    rel.matrix <- cached.fun(roc$boot.tpr, roc$boot.fpr, steps)
  } else {
    rel.matrix <- uncached.fun(roc$predictions, roc$true.classes, roc$n.boot, steps)
  }
  
  rm(roc)
  estimate <- data.frame(as.numeric(estimate))
  conf.area <- t(apply(rel.matrix, 2, quantile, alpha.levels))
  conf.area <- as.data.frame(conf.area)
  conf.area <- cbind(estimate, conf.area)
  conf.area <- cbind(data.frame(1 - seq(0, 1, by = (1 / steps))), conf.area)
  names(conf.area) <- frame.names
  class(conf.area) <- c("fbroc.conf", class(conf.area))
  return(conf.area)
}

#' Process bootstrapped TPR/FPR at thresholds matrix into TPR at FPR matrix
#' 
#' Usually \code{fbroc} calculates the TPR and FPR at the thresholds of the ROC curve.
#' per bootstrap replicate. This can be calculated quickly, but is often not convenient
#' to work with. Therefore \code{boot.tpr.at.fpr} instead gets the TPR along a sequence
#' of different values for the FPR.
#' @param roc Object of class \code{fbroc.roc}.
#' @param steps Number of discrete steps for the FPR at which the TPR is 
#' calculated. TPR confidence intervals are given for all FPRs in 
#' \code{seq(0, 1, by = (1 / steps))}. Defaults to \code{n.neg}, thus covering all possible values.
#' @return Matrix containing the TPR bootstrap replicates for the discrete
#' FPR steps.
#' @export
#' @seealso \code{\link{boot.roc}}
boot.tpr.at.fpr <- function(roc, steps = roc$n.neg) {
  steps = as.integer(steps)
  if (roc$use.cache) {
    rel.matrix <- tpr_at_fpr_cached(roc$boot.tpr, roc$boot.fpr, steps)
  } else {
    rel.matrix <- tpr_at_fpr_uncached(roc$predictions,
                                      as.integer(roc$true.classes),
                                      roc$n.boot,
                                      steps)
  }
  FPR.VEC = round(1 - seq(0, 1, by = (1 / steps)), 3)
  colnames(rel.matrix) <- paste("TPR.AT.FPR.", FPR.VEC, sep = "")
  return(rel.matrix)
}
