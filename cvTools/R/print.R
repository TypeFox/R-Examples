# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method print cvFolds
print.cvFolds <- function(x, ...) {
    # print general information
    if(x$n == x$K) {
        cvText <- "Leave-one-out CV"
    } else {
        cvText <- sprintf("%d-fold CV", x$K)
        if(x$R > 1) {
            cvText <- paste("Repeated", cvText, "with", x$R, "replications")
        }
    }
    cat(paste("\n", cvText, ":", sep=""))
    # print information on folds (add space between folds and subsets)
    subsets <- x$subsets
    if(x$R == 1) {
        cn <- "Index"
        nblanks <- 2
    } else {
        cn <- seq_len(x$R)
        nblanks <- 3
    }
    nblanks <- max(nchar(as.character(subsets[, 1])) - nchar(cn), 0) + nblanks
    cn[1] <- paste(c(rep.int(" ", nblanks), cn[1]), collapse="")
    dimnames(subsets) <- list(Fold=x$which, cn)
    print(subsets, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print cv
#' @S3method print summary.cv
print.cv <- print.summary.cv <- function(x, ...) {
    # print cross-validation results
    if(x$n == x$K) {
        cvText <- "Leave-one-out CV results:\n"
    } else cvText <- sprintf("%d-fold CV results:\n", x$K)
    cat(cvText)
    print(x$cv, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print cvSelect
#' @S3method print summary.cvSelect
print.cvSelect <- print.summary.cvSelect <- function(x, best = TRUE, ...) {
    # print cross-validation results
    if(length(K <- unique(x$K)) == 1) {
        if(x$n == K) {
            cat("\nLeave-one-out CV results:\n")
        } else cat(sprintf("\n%d-fold CV results:\n", K))
    } else cat("\nCV results:\n")
    print(x$cv, ...)
    # print optimal model if requested
    if(isTRUE(best)) {
        cat("\nBest model:\n")
        best <- x$best
        bestFit <- x$cv[best, "Fit"]
        if(is.factor(bestFit)) bestFit <- as.character(bestFit)
        names(bestFit) <- names(best)
        print(bestFit, ...)
    }
    # return object invisibly
    invisible(x)
}

#' @S3method print cvTuning
#' @S3method print summary.cvTuning
print.cvTuning <- print.summary.cvTuning <- function(x, best = TRUE, ...) {
    # print cross-validation results
    if(x$n == x$K) {
        cat("\nLeave-one-out CV results:\n")
    } else cat(sprintf("\n%d-fold CV results:\n", x$K))
    print(cbind(x$tuning, x$cv[, -1, drop=FALSE]), ...)
    # print optimal value for tuning parameters if requested
    if(isTRUE(best)) {
        if(ncol(x$tuning) == 1) {
            cat("\nOptimal tuning parameter:\n")
        } else cat("\nOptimal tuning parameters:\n")
        best <- x$best
        optimalTuning <- x$tuning[best, , drop=FALSE]
        rownames(optimalTuning) <- names(best)
        print(optimalTuning, ...)
    }
    # return object invisibly
    invisible(x)
}
