#' Resampling schemes
#' 
#' Performance evaluation and parameter tuning use resampling methods to
#' estimate the performance of models. These are defined by resampling
#' schemes, which are data frames where each column corresponds to a
#' division of the data set into mutually exclusive training and test sets.
#' Repeated hold out and cross-validation are two methods to create such
#' schemes.
#'
#' Note that when setting up analyzes, the user should not call
#' \code{resample_holdout} or \code{resample_crossvalidation} directly, as
#' \code{resample} performs additional necessary processing of the scheme.
#'
#' Resampling scheme can be visualized in a human digestible form with the
#' \code{\link[=image.resample]{image}} function.
#' 
#' Functions for generating custom resampling schemes should be implemented as
#' follows and then called by \code{resample("myMethod", ...)}:
#'
#' \code{resample_myMethod <- function(y, ..., subset)}
#' \describe{
#'     \item{\code{y}}{Response vector.}
#'     \item{\code{...}}{Method specific attributes.}
#'     \item{\code{subset}}{Indexes of observations to be excluded for the
#'         resampling.}
#' }
#' The function should return a list of the following elements:
#' \describe{
#'     \item{\code{folds}}{A data frame with the folds of the scheme that
#'         conforms to the description in the 'Value' section below.}
#'     \item{\code{parameter}}{A list with the parameters necessary to generate
#'         such a resampling scheme. These are needed when creating subschemes
#'         needed for parameter tuning, see \code{\link{subresample}}.}
#' }
#' 
#' @param method The resampling method to use, e.g. \code{"holdout"} or
#'   \code{"crossvalidation"}.
#' @param y Observations to be divided.
#' @param ... Sent to the method specific function, e.g.
#'   \code{"resample_holdout"}.
#' @param nfold Number of folds.
#' @param balanced Whether the sets should be balanced or not, i.e. if the
#'   class ratio over the sets should be kept constant (as far as possible).
#' @param subset Which objects in \code{y} that are to be divided and which
#'   that are not to be part of neither set.
#'   If \code{subset} is a resampling scheme, a list of inner
#'   cross-validation schemes will be returned.
#' @return A data frame defining a resampling scheme. \code{TRUE} or a positive integer
#'   codes for training set and \code{FALSE} or \code{0} codes for test set.
#'   Positive integers > 1 code for multiple copies of an observation in the
#'   training set. \code{NA} codes for neither training nor test set and is
#'   used to exclude observations from the analysis altogether.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{subresample}},
#'   \code{\link{image.resample}}, \code{\link{index_fit}}
#' @export
resample <- function(method, y, ..., subset=TRUE){
    stopifnot(is.character(method))
    if(inherits(y, "Surv")) y <- factor(dichotomize(y))

    method_fun <- match.fun(sprintf("resample_%s", method))
    res <- method_fun(y=y, ..., subset=subset)
    stopifnot(inherits(res$fold_set, "data.frame"))
    res$fold_set[T] <- Map(function(f, n){
        structure(f, class=c(method, "fold", class(f)),
                  parameter=res$parameter, fold.name=n)
    }, res$fold_set, names(res$fold_set))

    class(res$fold_set) <- c(method, "resample", class(res$fold_set))
    if(is.null(colnames(res$fold_set))){
        colnames(res$fold_set) <- sprintf("fold%i", 1:ncol(res$fold_set))
    }
    if(!is.null(names(y))) rownames(res$fold_set) <- names(y)

    # Confirm that the resampling method did not include any observations that
    # should be left out
    y_remove <- is.na(y) | !logical_subset(y, subset)
    if(any(y_remove)) res$fold_set[y_remove,] <- NA
    res$fold_set
}

#' Generate resampling subschemes
#'
#' A subscheme is a resampling scheme that only includes observations in the
#' training set of a fold. This function
#' automatically fetches the type and parameters of the prototype fold and use
#' them to generate the subscheme.
#'
#' @param fold A resampling scheme or fold to use to define the sub scheme(s).
#' @param y The observations used to create the resampling scheme. See
#'   \code{\link{resample}} for details.
#' @return A resampling scheme.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @examples
#' cv <- resample("holdout", y=1:12, test_fraction=1/4, nfold=3)
#' inner.cv <- subresample(cv, y=1:12)
#' @seealso \code{\link{emil}}, \code{\link{resample}}
#' @export
subresample <- function(fold, y){
    if(is.data.frame(fold)){
        # The user inputted a full scheme
        lapply(fold, subresample, y)
    } else {
        do.call(resample,
            c(list(method = class(fold)[1], y = y),
              attr(fold, "parameter"),
              list(subset = as.logical(fold))))
    }
}


#' Convert a fold to row indexes of fittdng or test set
#'
#' @param fold A fold of a resampling scheme.
#' @param allow_oversample Whether or not to allow individual observation to
#'   exist in multiple copies in the training set. This is typically not the
#'   case, but can be used when a class is underrepresented in the data set.
#' @return An integer vector of row indexes.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{resample}}
#' @export
index_fit <- function(fold, allow_oversample=TRUE){
    stopifnot(inherits(fold, "fold"))
    if(allow_oversample){
        rep(seq_along(fold), na_fill(fold, 0))
    } else {
        which(!fold %in% c(0, NA))
    }
}
#' @rdname index_fit
#' @aliases index_test
#' @export
index_test <- function(fold){
    stopifnot(inherits(fold, "fold"))
    which(fold %in% 0)
}


#' @param test_fraction Fraction of objects to hold out (0 < test_fraction < 1).
#' @examples
#' resample("holdout", 1:50, test_fraction=1/3)
#' resample("holdout", factor(runif(60) >= .5))
#' @rdname resample
#' @export
resample_holdout <- function(y, test_fraction=.5, nfold=5, balanced=is.factor(y), subset){
    res <- resample_bootstrap(y, nfold=nfold, fit_fraction=1-test_fraction, replace=FALSE,
                              balanced=balanced, subset=subset)
    res$parameter <- list(test_fraction = test_fraction, nfold = nfold,
                          balanced = balanced)
    res
}

#' @param nrepeat Number of fold sets to generate.
#' @examples
#' y <- factor(runif(60) >= .5)
#' cv <- resample("crossvalidation", y)
#' image(cv, main="Cross-validation scheme")
#' @rdname resample
#' @export
resample_crossvalidation <- function(y, nfold=5, nrepeat=5, balanced=is.factor(y), subset){
    n <- length(y)
    stopifnot(n >= nfold)
    subset <- positive_integer_subset(y, subset)

    fold_set <- as.data.frame(replicate(nrepeat, {
        idx <- if(!balanced){
            sample(subset)
        } else {
            levs <- if(is.factor(y)) levels(y) else unique(y)
            unlist(lapply(levs[order(table(y[subset]))], function(lev){
                sample(intersect(which(y == lev), subset))
            }))
        }
        idx <- matrix(c(idx, rep(NA, ceiling(length(idx)/nfold)*nfold-length(idx))),
                      ncol=nfold, byrow=TRUE)
        apply(idx, 2, function(i) !seq_len(n) %in% i)
    }))
    names(fold_set) <- sprintf("rep%ifold%i", rep(1:nrepeat, each=nfold),
                               rep(1:nfold, nrepeat))
    list(fold_set = fold_set,
         parameter = list(nfold=nfold, nrepeat=nrepeat, balanced=balanced))
}

#' @param fit_fraction The size of the training set relative to the entire data
#'   set.
#' @param replace Whether to sample with replacement.
#' @rdname resample
#' @export
resample_bootstrap <- function(y, nfold=10, fit_fraction = if(replace) 1 else .632,
                               replace=TRUE, balanced = is.factor(y), subset){
    stopifnot(fit_fraction <= 1 || replace)
    subset <- positive_integer_subset(y, subset)
    n <- length(subset)
    nsample <- round(fit_fraction*n)
    if(balanced){
        stopifnot(is.factor(y))
        class_n <- table(y[subset])
        class_idx <- split(subset, y[subset])
        class_sample <- fit_fraction * class_n

        class_min <- floor(class_sample)
        class_count <- rep(0, length(class_sample))
        class_rest <- class_sample-class_min
        npick <- nsample - sum(class_min)
        fold_sample <- structure(rep(list(class_min), nfold),
                                 names = sprintf("fold%i", 1:nfold))
        for(i in 1:nfold){
            class_count <- class_count + class_rest
            pick <- tail(order(class_count), npick)
            fold_sample[[i]][pick] <- fold_sample[[i]][pick] + 1
            class_count[pick] <- class_count[pick] - 1
        }
        fold_set <- lapply(fold_sample, function(fs){
            Map(function(idx, n){
                sample(idx, n, replace=replace)
            }, class_idx, fs) %>%
                unlist %>% factor(seq_along(y)) %>% table %>% unclass
        }) %>%
            as.data.frame %>%
            return
    } else {
        fold_set <- replicate(nfold, {
            sample(subset, nsample, replace=replace) %>%
                unlist %>% factor(seq_along(y)) %>% table %>% unclass
        }, simplify = FALSE) %>%
            as.data.frame %>%
            structure(names=sprintf("fold%i", 1:nfold)) %>%
            return
    }
    list(fold_set = fold_set,
         parameter = list(nfold=nfold, fit_fraction=fit_fraction,
                          replace=replace, balanced=balanced))
}

#' Visualize resampling scheme
#'
#' Class specific extension to \code{\link{image}}.
#' 
#' @method image resample
#' @param x Resampling scheme, as returned by \code{\link{resample}}.
#' @param col Color palette matching the values of \code{x}.
#'   Can also be the response vector used to create the
#'   scheme for automatic coloring.
#' @param ... Sent to \code{\link{plot}}.
#' @return Nothing, produces a plot.
#' @examples
#' y <- gl(2, 30)
#' image(resample("crossvalidation", y, nfold=3, nrepeat=8), col=y)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{resample}}
#' @export
image.resample <- function(x, col, ...){
    x <- as.matrix(x)
    if(missing(col)) col <- gl(1, nrow(x))
    if(inherits(col, "Surv")) col <- dichotomize(col, to_factor=TRUE)
    if(is.factor(col)){
        y <- col
        if(length(y) != nrow(x))
            stop("Color vector does not match resampling scheme.")
        nice_require("RColorBrewer")
        if(length(levels(y)) > 12)
            warning("Too few colors to assign unique ones to each class.")
        col <- rep(RColorBrewer::brewer.pal(12, "Set3"),
                   ceiling(length(levels(y))/12))[seq_along(levels(y))]
        col <- apply(col2rgb(col), 2, function(cl){
              apply(cl %o% seq(.7, 1, length.out=max(x, na.rm=TRUE)+1), 2,
                    function(cc) do.call(rgb, as.list(cc/255)))
        })
        mat <- matrix(col[cbind(as.vector(x)+1, as.integer(y))], nrow(x), ncol(x))
    } else {
        mat <- 1 + x
        mat <- matrix(col[mat], nrow(mat))
    }
    mat[is.na(x)] <- "transparent"
    par(xaxs="i", yaxs="i")
    plot(c(.5, ncol(x)+.5), c(.5, nrow(x)+.5), type="n",
         axes=FALSE, xlab="Folds", ylab="Observations", ...)
    rasterImage(mat, .5, .5, ncol(x)+.5, nrow(x)+.5, interpolate=FALSE)
    nice_axis(1)
    ticks <- pretty(pretty(par("usr")[3:4]))
    nice_axis(2, at=nrow(x)-ticks+1, labels=ticks, las=1)
    nice_box()
}


#' @method image crossvalidation
#' @rdname image.resample
#' @export
image.crossvalidation <- function(x, col, ...){
    image.resample(x, col, ...)
    parameter <- attr(x[[1]], "parameter")
    if(parameter$nrepeat > 1){
        l <- 1:(parameter$nrepeat-1)*parameter$nfold + .5
        segments(l, par("usr")[3], l, par("usr")[4])
    }
}

#' Extract parts of a resampling scheme
#' 
#' @param x Resampling scheme or fold, as created by \code{\link{resample}}.
#' @param ... Sent to a class specific extraction function, see \link{Extract}.
#' @examples
#' cv <- resample("crossvalidation", iris$Species)
#' cv[1:4,]
#' cv[, 2:7]
#' cv[1:4, 2:7]
#' cv[,1]
#' cv[1:4, 1, drop=FALSE]
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
#' @export
`[.resample` <- function(x, ...){
    original_class <- class(x)
    class(x) <- setdiff(original_class, "resample")
    x <- x[...]
    if(!inherits(x, "fold")){
        class(x) <- original_class
    }
    x
}
#' @noRd
#' @export
`[.fold` <- function(x, ...){
    attribute <- attributes(x)
    class(x) <- setdiff(class(x), "fold")
    x <- x[...]
    attributes(x) <- attribute
    x
}

