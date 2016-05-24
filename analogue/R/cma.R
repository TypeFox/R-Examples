###########################################################################
##                                                                       ##
## cma           - extracts and formats close modern analogues           ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object for method dispatch. Only class 'analog'.      ##
## cutoff        - numeric. Critical value determining level above which ##
##                 samples are defined as close modern analogues         ##
##                                                                       ##
###########################################################################
cma <- function(object, ...) UseMethod("cma")

cma.default <- function(object, ...)
  {
    stop("No default method for \"cma\"")
  }

cma.analog <- function(object, cutoff, prob = c(0.01, 0.025, 0.05), ...)
  {
    if (!inherits(object, "analog"))
      stop("Use only with \"analog\" objects")
    if(missing(cutoff)) {
      if(is.null(object$train))
        stop("If 'cutoff' is not provided, 'object' must contain\ncomponent \"train\"")
      else
        cutoff <- quantile(dissim(object), probs = 0.025)
    } else {
      if (!is.numeric(cutoff))
        stop("Argument \"cutoff\" must be numeric")
    }
    #if(!any(apply(object$analogs, 2, function(x) any(x <= cutoff))))
    #  stop(paste("No analogues as close or closer than \"cutoff = ",
    #             cutoff, "\":\n\tChoose a more suitable value", sep = ""))
    n.samp <- ncol(object$analogs)
    nams <- colnames(object$analogs)
    ## don't want apply() as that may simplify if all samples have
    ## same number of analogues - must return a list
    nc <- ncol(object$analogs)
    close <- vector("list", nc)
    for(i in seq_len(nc)) {
        close[[i]] <- sortByCutoff(object$analogs[, i], cutoff)
    }
    names(close) <- colnames(object$analogs)
    if(length(close) == 0)
      close <- vector(mode = "list", length = length(nams))
    each.analogs <- sapply(close, length, USE.NAMES = FALSE)
    names(each.analogs) <- names(close) <- nams
    .call <- match.call()
    .call[[1]] <- as.name("cma")
    structure(list(close = close,
                   call = .call, cutoff = cutoff,
                   quant = quantile(dissim(object), probs = prob),
                   prob = prob,
                   method = object$method,
                   n.analogs = each.analogs),
              class = "cma")
  }

## First attempt at this method - we want k to select the k closest analogues
## but also allow cutoff for later when mat will work with a threshold
`cma.mat` <- function(object, k, cutoff, prob = c(0.01, 0.025, 0.05), ...) {
    if (!inherits(object, "mat"))
        stop("Use only with \"mat\" objects")
    n.samp <- ncol(object$Dij)
    nams <- colnames(object$Dij)
    K <- !missing(k)
    CUT <- !missing(cutoff)
    if(K && CUT)
        stop("Only one of \"k\" and \"cutoff\" may be used, not both.")
    if(!K && !CUT) {
        k <- getK(object)
        cutoff <- NULL
        K <- TRUE
    }
    if(K) {
        close <- vector(mode = "list", length = n.samp)
        ks <- seq_len(k)
        for(i in seq_along(close)) {
            close[[i]] <- sortByK(object$Dij[, i], ks)
        }
        each.analogs <- sapply(close, length, USE.NAMES = FALSE)
        names(each.analogs) <- names(close) <- nams
    } else {
        ## don't want apply as that may simplify if all samples have
        ## same number of analogues - must return a list
        nc <- ncol(object$Dij)
        close <- vector("list", nc)
        for(i in seq_len(nc)) {
            close[[i]] <-  sortByCutoff(object$Dij[,i], cutoff)
        }
        names(close) <- colnames(object$Dij)
        each.analogs <- sapply(close, length, USE.NAMES = FALSE)
        k <- NULL
        names(each.analogs) <- names(close) <- nams
    }
    if(length(close) == 0) {
        close <- vector(mode = "list", length = length(nams))
        names(each.analogs) <- names(close) <- nams
    }
    .call <- match.call()
    .call[[1]] <- as.name("cma")
    structure(list(close = close,
                   call = .call, cutoff = cutoff, k = k,
                   quant = quantile(dissim(object), probs = prob,
                   na.rm = TRUE),
                   prob = prob,
                   method = object$method,
                   n.analogs = each.analogs),
              class = "cma")
}
## First attempt at this method - we want k to select the k closest analogues
## but also allow cutoff for later when mat will work with a threshold
`cma.predict.mat` <- function(object, k, cutoff, prob = c(0.01, 0.025, 0.05),
                              ...) {
    if (!inherits(object, "predict.mat"))
        stop("Use only with \"predict.mat\" objects")
    n.samp <- ncol(object$Dij)
    nams <- colnames(object$Dij)
    K <- !missing(k)
    CUT <- !missing(cutoff)
    if(K && CUT)
        stop("Only one of \"k\" and \"cutoff\" may be used, not both.")
    if(!K && !CUT) {
        k <- getK(object)
        cutoff <- NULL
        K <- TRUE
    }
    if(K) {
        close <- vector(mode = "list", length = n.samp)
        ks <- seq_len(k)
        for(i in seq_along(close)) {
            close[[i]] <- sortByK(object$Dij[, i], ks)
        }
        each.analogs <- sapply(close, length, USE.NAMES = FALSE)
        names(each.analogs) <- names(close) <- nams
    } else {
        ## don't want apply as that may simplify if all samples have
        ## same number of analogues - must return a list
        nc <- ncol(object$Dij)
        close <- vector("list", nc)
        for(i in seq_len(nc)) {
            close[[i]] <-  sortByCutoff(object$Dij[,i], cutoff)
        }
        names(close) <- colnames(object$Dij)
        each.analogs <- sapply(close, length, USE.NAMES = FALSE)
        k <- NULL
        names(each.analogs) <- names(close) <- nams
    }
    if(length(close) == 0) {
        close <- vector(mode = "list", length = length(nams))
        names(each.analogs) <- names(close) <- nams
    }
    .call <- match.call()
    .call[[1]] <- as.name("cma")
    structure(list(close = close,
                   call = .call, cutoff = cutoff, k = k,
                   quant = NULL,
                   prob = prob,
                   method = object$method,
                   n.analogs = each.analogs),
              class = "cma")
}

print.cma <- function(x, digits = min(3, getOption("digits") - 4), ...) {
    method <- x$method
    .call <- deparse(x$call)
    cat("\n")
    writeLines(strwrap("Close modern analogues of fossil samples",
                       prefix = "\t"))
    cat(paste("\nCall:", .call, "\n"))
    cat(paste("\nDissimilarity:", method, "\n\n"))
    if(is.null(x$k)) {
        cat("     k: Not supplied\n\n")
    } else {
        cat(paste("     k:", x$k[1], "\n"))
    }
    if(is.null(x$cutoff)) {
        cat("Cutoff: Not supplied\n\n")
    } else {
        cat(paste("Cutoff:", round(x$cutoff, digits), "\n\n"))
    }
    writeLines(strwrap("Number of analogues per fossil sample:",
                       prefix = "\t"))
    cat("\n")
    print(x$n.analogs, digits = digits)
    cat("\n")
    invisible(x)
}

## two simple functions that now get used more often in the methods above
## so made full functions (rather than in-line) but not exported so
## don't need to be documented
sortByCutoff <- function(x, cutoff) {
    x <- sort(x)
    x <- x[x <= cutoff]
}

sortByK <- function(x, ks) {
    x <- sort(x)
    x[ks]
}
