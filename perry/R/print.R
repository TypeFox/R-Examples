# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @S3method print cvFolds
print.cvFolds <- function(x, ...) {
  # print general information
  cvText <- getPrefix(x)
  if(x$R > 1) cvText <- paste(cvText, "with", x$R, "replications")
  cat(paste("\n", cvText, ":", sep=""))
  # print information on folds (add space between folds and subsets)
  subsets <- x$subsets
  if(x$R == 1) {
    cn <- if(hasComponent(x, "grouping")) "Group index" else "Index"
    nblanks <- 2
  } else {
    cn <- as.character(seq_len(x$R))
    nblanks <- 3
  }
  nblanks <- max(nchar(as.character(subsets[, 1]))-nchar(cn[1]), 0) + nblanks
  cn[1] <- paste(c(rep.int(" ", nblanks), cn[1]), collapse="")
  dimnames(subsets) <- list(Fold=x$which, cn)
  print(subsets, ...)
  # return object invisibly
  invisible(x)
}

#' @S3method print randomSplits
print.randomSplits <- function(x, ...) {
  # print general information
  if(x$R == 1) {
    cat("\nRandom split\n")
    cn <- if(hasComponent(x, "grouping")) "Group index" else "Index"
  } else {
    cat(sprintf("\n%d random splits\n", x$R))
    cn <- as.character(seq_len(x$R))
  }
  prefix <- if(hasComponent(x, "grouping")) "Groups" else "Observations"
  cat(sprintf("%s in test data:\n", prefix))
  # print indices of items in test data
  subsets <- x$subsets
  colnames(subsets) <- cn
  print(subsets, ...)
  # return object invisibly
  invisible(x)
}

#' @S3method print bootSamples
print.bootSamples <- function(x, ...) {
  # print general information
  if(x$R == 1) {
    postfix <- "sample"
    cn <- if(hasComponent(x, "grouping")) "Group index" else "Index"
  } else {
    postfix <- "samples"
    cn <- as.character(seq_len(x$R))
  }
  prefix <- if(hasComponent(x, "grouping")) "Groups" else "Observations"
  cat(sprintf("\n%s in the bootstrap %s:\n", prefix, postfix))
  # print indices of items in test data
  samples <- x$samples
  colnames(samples) <- cn
  print(samples, ...)
  # return object invisibly
  invisible(x)
}

#' @S3method print perry
#' @S3method print summary.perry
print.perry <- print.summary.perry <- function(x, ...) {
  # print cross-validation results
  cat(getPrefix(x$splits), "results:\n")
  print(x$pe, ...)
  # return object invisibly
  invisible(x)
}

#' @S3method print perrySelect
#' @S3method print summary.perrySelect
print.perrySelect <- print.summary.perrySelect <- function(x, best = TRUE, ...) {
  # print cross-validation results
  cat("\n")
  cat(getPrefix(x$splits), "results:\n")
  print(x$pe, ...)
  # print optimal model if requested
  if(isTRUE(best)) {
    cat("\nBest model:\n")
    best <- x$best
    bestFit <- x$pe[best, "Fit"]
    if(is.factor(bestFit)) bestFit <- as.character(bestFit)
    names(bestFit) <- names(best)
    print(bestFit, ...)
  }
  # return object invisibly
  invisible(x)
}

#' @S3method print perryTuning
#' @S3method print summary.perryTuning
print.perryTuning <- print.summary.perryTuning <- function(x, best = TRUE, 
                                                           final = TRUE, ...) {
  # print cross-validation results
  cat("\n")
  cat(getPrefix(x$splits), "results:\n")
  print(cbind(x$tuning, x$pe[, -1, drop=FALSE]), ...)
  # print optimal value for tuning parameters if requested
  if(isTRUE(best)) {
    if(ncol(x$tuning) == 1) cat("\nOptimal tuning parameter:\n")
    else cat("\nOptimal tuning parameters:\n")
    best <- x$best
    optimalTuning <- x$tuning[best, , drop=FALSE]
    rownames(optimalTuning) <- names(best)
    print(optimalTuning, ...)
  }
  # print final model
  if(isTRUE(final) && !is.null(finalModel <- x$finalModel)) {
    cat("\nFinal model:\n")
    print(finalModel, ...)
  }
  # return object invisibly
  invisible(x)
}


## get prefix for print methods

getPrefix <- function(x) UseMethod("getPrefix")

getPrefix.cvFolds <- function(x) {
  if(x$n == x$K) "Leave-one-out CV"
  else sprintf("%d-fold CV", x$K)
}

getPrefix.randomSplits <- function(x) "Random splitting"

getPrefix.bootSamples <- function(x) {
  if(x$type == "out-of-bag") "Out-of-bag bootstrap"
  else sprintf("%s bootstrap", x$type)
}

getPrefix.default <- function(x) "Prediction error"
