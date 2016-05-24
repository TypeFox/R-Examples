#' Bootstrap paired ROC curves
#'
#' Given two numerical predictors for the same outcome on the same set of samples, this functions
#' enables the bootstrapping of the paired ROC curves of the two prediction models. While bootstrapping
#' the same set of samples are used for both curves in each iteration, preserving the correlation
#' between the two models.  
#'
#' @inheritParams boot.roc
#' @param pred1 Numerical predictions for the first classifier.
#' @param pred2 Numerical predictions for the second classifier.
#' @return A list of class \code{fbroc.paired.roc}, containing the elements:
#' \item{prediction1}{Input predictions for first model.}
#' \item{prediction2}{Input predictions for second model.}
#' \item{true.class}{Input classes.}
#' \item{n.thresholds1}{Number of thresholds of the first predictor.}
#' \item{n.thresholds2}{Number of thresholds of the second predictor.}
#' \item{n.boot}{Number of bootstrap replicates.}
#' \item{use.cache}{Indicates if cache is used for this ROC object.}
#' \item{tie.strategy}{Used setting how to handle ties in predictors.}
#' \item{n.pos}{Number of positive observations.}
#' \item{n.neg}{Number of negative observations.}
#' \item{roc1}{A data.frame containing the thresholds of the first ROC curve and the TPR and FPR at these
#' thresholds.}
#' \item{roc2}{A data.frame containing the thresholds of the second ROC curve and the TPR and FPR at these
#' thresholds.}
#' \item{auc1}{The AUC of the first ROC curve.}
#' \item{auc2}{The AUC of the second ROC curve.}
#' \item{boot.tpr1}{If the cache is enabled, a matrix containing the bootstrapped TPR at the thresholds
#' for the first predictor.}
#' \item{boot.fpr1}{If the cache is enabled, a matrix containing the bootstrapped FPR at the thresholds
#' for the first predictor.}
#' \item{boot.tpr2}{If the cache is enabled, a matrix containing the bootstrapped TPR at the thresholds
#' for the second predictor.}
#' \item{boot.fpr2}{If the cache is enabled, a matrix containing the bootstrapped FPR at the thresholds
#' for the second predictor.}
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
#' data(roc.examples)
#' # Do not use cache
#' example <- boot.paired.roc(roc.examples$Cont.Pred, roc.examples$Cont.Pred.Outlier,
#'                           roc.examples$True.Class, n.boot = 500)
#' perf(example, "auc") # estimate difference in auc
#' perf(example, "tpr", fpr = 0.5) # estimate difference in TPR at a FPR of 50%
#' plot(example) # show plot
#' # Cached mode
#' example <- boot.paired.roc(roc.examples$Cont.Pred, roc.examples$Cont.Pred.Outlier,
#'                           roc.examples$True.Class, n.boot = 1000, use.cache = TRUE)
#' conf(example, conf.for = "tpr", steps = 10) # get confidence regions for TPR at FPR
#' conf(example, conf.for = "fpr", steps = 10) # get confidence regions for FPR at TPR
#' perf(example, "fpr", tpr = 0.9) # estimate difference in FPR at a TPR of 90%                     
#' @seealso \code{\link{boot.roc}}, 
#'          \code{\link{plot.fbroc.paired.roc}},
#'          \code{\link{perf.fbroc.paired.roc}} 
#' 
#' @export
boot.paired.roc <- function(pred1, pred2, true.class, stratify = TRUE, n.boot = 1000,
                            use.cache = FALSE, tie.strategy = NULL) {
  # validate input
  if ((length(pred1) != length(true.class)))
    stop("Predictions and true classes need to have the same length")
  if ((length(pred2) != length(true.class)))
    stop("Predictions and true classes need to have the same length")
  if (class(pred1) == "integer") {
    pred1 <- as.numeric(pred1)
  }
  if ((class(pred1) != "numeric"))
    stop("Predictions must be numeric")
  
  if (class(pred2) == "integer") {
    pred2 <- as.numeric(pred2)
  }
  if ((class(pred2) != "numeric"))
    stop("Predictions must be numeric")
  
  if ((class(true.class) != "logical"))
    stop("Classes must be logical")
  if ((class(stratify) != "logical"))
    stop("Classes must be logical")
  
  index.na <- is.na(pred1) | is.na(pred2) | is.na(true.class)
  if (any(index.na)) {
    n <- sum(index.na)
    warning.msg <- 
      paste(n, "observations had to be removed due to missing values")
    warning(warning.msg)
    true.class <- true.class[!index.na]
    pred1 <- pred1[!index.na]
    pred2 <- pred2[!index.na]
  }
  if (is.null(tie.strategy)) {
    if ((length(pred1) == length(unique(pred1))) &
        (length(pred2) == length(unique(pred2)))) 
      tie.strategy <- 1 
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
  original.rocs <- paired_roc_analysis(pred1, pred2, true.int)
  
  roc1 = c.list.to.roc(original.rocs[[1]])
  roc2 = c.list.to.roc(original.rocs[[2]])
  
  auc1 = original.rocs[[1]][[4]]
  auc2 = original.rocs[[2]][[4]]
  if (use.cache) {
    
    booted.roc <- tpr_fpr_boot_paired(pred1, pred2, true.int, n.boot)
    boot.tpr1 <- booted.roc[[1]]
    boot.fpr1 <- booted.roc[[2]]
    boot.tpr2 <- booted.roc[[3]]
    boot.fpr2 <- booted.roc[[4]]
    rm(booted.roc)
  } else {
    boot.tpr1 <- NULL
    boot.fpr1 <- NULL
    boot.tpr2 <- NULL
    boot.fpr2 <- NULL
  }
  
  output <- list(predictions1 = pred1,
                 predictions2 = pred2,
                 true.classes = true.class,
                 n.thresholds1 = nrow(roc1),
                 n.thresholds2 = nrow(roc2),
                 n.boot = as.integer(n.boot),
                 use.cache = use.cache,
                 tie.strategy = tie.strategy,
                 n.pos = sum(true.class),
                 n.neg = sum(!true.class),
                 roc1 = roc1,
                 roc2 = roc2,
                 auc1 = auc1,
                 auc2 = auc2,
                 #auc = auc,
                 boot.tpr1 = boot.tpr1,
                 boot.fpr1 = boot.fpr1,
                 boot.tpr2 = boot.tpr2,
                 boot.fpr2 = boot.fpr2)
  class(output) <- append(class(output), "fbroc.paired.roc")
  return(output)
}

#' Generates confidence intervals for the difference in TPR between two predictors for a range of FPRs or vice versa
#' 
#' Calculates confidence intervals for the difference in TPR at different FPR values or vice versa. The stepsize
#' at which the TPR or FPR is calculated can be set as needed.
#' 
#' @param roc An object of class \code{fbroc.paired.roc}.
#' @inheritParams conf.fbroc.roc
#' @return A data.frame containing either discrete TPR steps and estimates and confidence bounds for
#' the difference FPR or vice versa, depending upon \code{conf.for}.
#' @section Details:
#' This function only gives estimates and confidence for the difference in the requested rate
#' (either True Positive Rate or False Positive Rate) between the first and the second classifier.
#' The values given are positive iff the first classifier has a higher rate. To get confidence regions
#' for either of the two underlying ROC curves use \code{conf} on the result of \code{extract.roc}.
#' @examples 
#' data(roc.examples)
#' example <- boot.paired.roc(roc.examples$Cont.Pred, roc.examples$Cont.Pred.Outlier,
#'                            roc.examples$True.Class, n.boot = 100)
#' conf(example, conf.for = "tpr", steps = 10) # get confidence regions for Delta TPR at FPR
#' conf(example, conf.for = "fpr", steps = 10) # get confidence regions for Delta FPR at TPR
#' @seealso \code{\link{boot.paired.roc}}, \code{\link{conf.fbroc.roc}},\code{\link{extract.roc}}
#' @export
conf.fbroc.paired.roc <- function(roc, conf.level = 0.95, conf.for = "TPR", steps = 250, ...) {
  conf.for <- toupper(conf.for)
  if (!(conf.for %in% c("TPR", "FPR"))) stop("Invalid rate given for confidence region")
  
  alpha <- 0.5*(1 - conf.level)
  alpha.levels <- c(alpha, 1 - alpha) 
  steps = as.integer(steps)
  
  if (conf.for == "TPR") {
    cached.fun <- tpr_at_fpr_delta_cached
    uncached.fun <- tpr_at_fpr_delta_uncached
    frame.names <- c("FPR", "Delta.TPR", "Lower.Delta.TPR", "Upper.Delta.TPR")
  } else {
    cached.fun <- fpr_at_tpr_delta_cached
    uncached.fun <- fpr_at_tpr_delta_uncached
    frame.names <- c("TPR", "Delta.FPR", "Lower.Delta.FPR", "Upper.Delta.FPR")
  }
  
  estimate <- cached.fun(matrix(roc$roc1$TPR, nrow = 1), 
                         matrix(roc$roc1$FPR, nrow = 1), 
                         matrix(roc$roc2$TPR, nrow = 1), 
                         matrix(roc$roc2$FPR, nrow = 1),
                         steps)
  
  if (roc$use.cache) {
    rel.matrix <- cached.fun(roc$boot.tpr1,
                             roc$boot.fpr1,
                             roc$boot.tpr2,
                             roc$boot.fpr2,
                             steps)
  } else {
    rel.matrix <- uncached.fun(roc$predictions1,
                               roc$predictions2,
                               as.integer(roc$true.classes),
                               roc$n.boot,
                               steps)
  }
  
  rm(roc)
  estimate <- data.frame(as.numeric(estimate))
  conf.area <- t(apply(rel.matrix, 2, quantile, alpha.levels))
  conf.area <- as.data.frame(conf.area)
  conf.area <- cbind(estimate, conf.area)
  conf.area <- cbind(data.frame(1 - seq(0, 1, by = (1 / steps))), conf.area)
  names(conf.area) <- frame.names
  class(conf.area) <- c("fbroc.conf.paired", class(conf.area))
  return(conf.area)
}


