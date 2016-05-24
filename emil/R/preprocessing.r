#' Data preprocessing
#' 
#' These functions are run in \code{\link{evaluate}} just prior to model
#' fitting, to extract fitting and test sets from the entire dataset and apply
#' transformations to pre-process the data (for handling missing values,
#' scaling, compression etc.).
#' They can also be used to adapt the form of the data to a specific
#' fitting function, e.g. \code{\link{pre_pamr}} that transposes the dataset
#' to make it compatible with the \code{pamr} classification method.
#' 
#' When supplied to \code{\link{evaluate}}, pre-processing functions can be
#' chained (i.e. executed sequentially) after an initating call to
#' \code{\link{pre_split}}.
#' This can either be done using the \code{\link[=chain]{pipe operator}} defined
#' in the \pkg{magrittr} package or by putting all pre-processing functions in a
#' regular list (see the examples).
#' 
#' Note that all transformations are defined based on the fitting data only
#' and then applied to both fitting set and test set. It is important to not let
#' the test data in any way be part of the model fitting, including the
#' preprocessing, to not risk information leakage and biased results!
#'
#' The imputation functions can also be used outside of
#' \code{\link{evaluate}} by not supplying a fold to
#' \code{\link{pre_split}}.
#' See the code of \code{\link{impute_median}} for an example.
#' 
#' @return A list with the following components
#' \describe{
#'     \item{\code{fit}}{Fitting set.}
#'     \item{\code{test}}{Test set.}
#'     \item{\code{feature_selection}}{Integer vector mapping the features of
#'           the training and test sets to the original data sets.}
#'     \item{\code{fold}}{The fold that was used to split the data.}
#' }
#'
#' @example examples/pre-process.r
#' @seealso \code{\link{pre_factor_to_logical}}, \code{\link{emil}},
#'   \code{\link{pre_impute_knn}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @name pre_process
{}

#' @param x Dataset.
#' @param y Response vector.
#' @param fold A logical or numeric vector with \code{TRUE} or positive numbers
#'   for fitting observations, \code{FALSE} or \code{0} for test
#'   observations, and \code{NA} for observations not to be included.
#' @rdname pre_process
#' @export
pre_split <- function(x, y, fold){
    if(missing(fold)){
        fold <- structure(rep(TRUE, nrow(x)),
                          class=c("fit_only", "fold", "logical"))
    }
    if(is.character(y) && length(y) == 1){
        y_col <- match(y, colnames(x))
        if(is.na(y_col))
            stop(paste("The is no column named", y, "in `x`."))
        y <- x[,y_col]
        x <- x[,-y_col]
    }
    structure(
        list(name = deparse(substitute(x)),
             fit = list(x = x[index_fit(fold),,drop=FALSE],
                        y = y[index_fit(fold)]),
             test = list(x = x[index_test(fold),,drop=FALSE],
                         y = y[index_test(fold)]),
             feature_selection = structure(1:ncol(x), names = colnames(x)),
             fold = fold),
        class=c("preprocessed_data", "list"))
}

#' Validate a pre-processed data set
#' 
#' While writing and debugging pre-processing functions this function can be
#' useful to confirm that the resulting data sets fulfills the necessary
#' requirements.
#' 
#' @param data Pre-processed data set.
#' @return Nothing, only throws an error or prints a completion message.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export 
validate_data <- function(data){
    stopifnot(inherits(data, "preprocessed_data"))
    stopifnot(identical(colnames(data$fit$x), colnames(data$test$x)))
    stopifnot(identical(colnames(data$fit$x), names(data$feature_selection)))
    cat("Data set is valid.")
}

#' Print method for pre-processed data
#' 
#' @method print preprocessed_data
#' @param x Pre-processed data, as produced by \code{\link{pre_split}}.
#' @param ... Ignored, kept for S3 consistency.
#' @return Nothing
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
print.preprocessed_data <- function(x, ...){
    feature <- table(factor(x$feature_selection, c(FALSE, TRUE)))
    cat("Pre-processed data set `", x$name, "` of ",
        ncol(x$fit$x), " features\n",
        nrow(x$fit$x), " observations for model fitting,\n",
        nrow(x$test$x), " observations for model evaluation.\n", sep="")
}

#' @param data Fitting and testing data sets, as returned by
#'   \code{\link{pre_split}}.
#' @param x_fun Function to apply to the descriptors of the datasets
#'   (e.g. \code{x}). This function will be applied independenly to the fitting
#'   and testing sets.
#' @param y_fun Function to be applied to the response of the training and test
#'   sets (independently).
#' @param ... Sent to internal methods, see the code of each function.
#' @rdname pre_process
#' @export
pre_convert <- function(data, x_fun, y_fun, ...){
    if(!missing(x_fun)){
        data$fit$x <- x_fun(data$fit$x, ...)
        data$test$x <- x_fun(data$test$x, ...)
    }
    if(!missing(y_fun)){
        data$fit$y <- y_fun(data$fit$y, ...)
        data$test$y <- y_fun(data$test$y, ...)
    }
    data
}

#' @rdname pre_process
#' @export
pre_transpose <- function(data){
    pre_convert(data, x_fun=t)
}

#' @param feature The features to be removed. Can be integer, logical or
#'   character.
#' @rdname pre_process
#' @export
pre_remove <- function(data, feature){
    if(!any(feature)){
        return(data)
    }
    if(is.logical(feature)){
        feature <- which(feature)
    } else if(is.character(feature)){
        feature <- match(feature, colnames(data$fit$x))
    }
    data$fit$x <- data$fit$x[, -feature, drop=FALSE]
    data$test$x <- data$test$x[, -feature, drop=FALSE]
    data$feature_selection <- data$feature_selection[-feature]
    stopifnot(ncol(data$fit$x) == length(data$feature_selection),
              ncol(data$test$x) == length(data$feature_selection))
    data
}

#' @param na.rm A logical value indicating whether \code{NA} values should be
#'   ignored.
#' @rdname pre_process
#' @export
pre_center <- function(data, y=FALSE, na.rm=TRUE){
    m <- colMeans(data$fit$x, na.rm=na.rm)
    data$fit$x <- sweep(data$fit$x, 2, m, "-")
    data$test$x <- sweep(data$test$x, 2, m, "-")
    if(y){
        my <- mean(y, na.rm=na.rm)
        if(is.na(my)) stop("Could not calculate the mean of `y`.")
        data$fit$y <- data$fit$y - my
        data$test$y <- data$test$y - my
    }
    data
}

#' @param center Whether to center the data before scaling.
#' @rdname pre_process
#' @export
pre_scale <- function(data, y=FALSE, na.rm=TRUE, center=TRUE){
    if(center) data <- pre_center(data, y=y)
    s <- apply(data$fit$x, 2, sd, na.rm=na.rm)
    data$fit$x <- sweep(data$fit$x, 2, s, "/")
    data$test$x <- sweep(data$test$x, 2, s, "/")
    if(y){
        sy <- sd(y)
        data$fit$y <- data$fit$y / sy
        data$test$y <- data$test$y / sy
    }
    data
}

#' @rdname pre_process
#' @export
pre_remove_constant <- function(data, na.rm=TRUE){
    if(is.data.frame(data$fit$x)){
        constant_feature <- vapply(data$fit$x, is_constant, logical(1), na.rm=na.rm)
    } else constant_feature <- apply(data$fit$x, 2, is_constant, na.rm=na.rm)
    if(anyNA(constant_feature)){
        if(na.rm) constant_feature[is.na(constant_feature)] <- TRUE
        else stop("Could not determine which features are constant.")
    } 
    pre_remove(data, constant_feature)
}

#' @param cutoff See \code{\link[caret]{findCorrelation}}.
#' @rdname pre_process
#' @export
pre_remove_correlated <- function(data, cutoff){
    if(missing(cutoff)) stop("`pre_remove_correlated` requires a cutoff.")
    nice_require("caret")
    pre_remove(data, caret::findCorrelation(cor(data$fit$x), cutoff = .75))
}

#' @param ncomponent Number of PCA components to use. Missing all components
#'   are used.
#' @param scale. Sent to \code{\link{prcomp}}.
#' @rdname pre_process
#' @export
pre_pca <- function(data, ncomponent, scale. = TRUE, ...){
    if(missing(ncomponent)){
        pca <- prcomp(data$fit$x, scale. = scale., ..., retx = TRUE)
        data$fit$x <- pca$x
    } else {
        pca <- prcomp(data$fit$x, scale. = scale., ..., retx = FALSE)
        pca$rotation <- pca$rotation[,1:ncomponent]
        data$fit$x <- predict(pca, data$fit$x)
    }
    data$test$x <- predict(pca, data$test$x)
    data
}

#' Convert factors to logical columns
#' 
#' Factors will be converted to one logical column per level (or one fewer if a
#' base level is specified).
#' 
#' @param data Pre-processed data set, as produced by \code{\link{pre_split}}.
#' @param feature Character vector with names of features to convert.
#'   Defaults to all factors in the data set.
#' @param base Sent to \code{\link{factor_to_logical}}. To specify different bases for
#'   different columns supply a vector or list with named elements.
#' @param drop Sent to \code{\link{factor_to_logical}}. To specify different bases for
#'   different columns supply a vector or list with named elements.
#' @examples
#' x <- mtcars[-1]
#' x <- transform(x,
#'     cyl = factor(cyl, ordered=TRUE),
#'     vs = factor(vs),
#'     gear = factor(gear)
#' )
#' y <- mtcars$mpg
#' cv <- resample("crossvalidation", y)
#' data <- pre_split(x, y, cv[[1]]) %>%
#'     pre_factor_to_logical(base = c(cyl="4", vs="0"), 
#'                           drop=c(cyl=FALSE, gear=FALSE))
#' data$fit$x
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pre_factor_to_logical <- function(data, feature, base=1L, drop=TRUE){
    stopifnot(is.data.frame(data$fit$x))
    convert <- function(x, base, base_missing, drop, drop_missing, name, warn=TRUE){
        withCallingHandlers({
            if(name %in% feature){
                factor_to_logical(x,
                    base = if(base_missing) 1L else base,
                    drop = if(drop_missing) TRUE else drop)
            } else x
        }, warning = function(w){
            if(warn){
                warning(sprintf("Column `%s` threw a warning when converted to logical: %s",
                                name, w$message))
            }
            invokeRestart("muffleWarning")
        }, error = function(e){
            stop(sprintf("Column `%s` could not be converted to logical: %s",
                         name, e$message))
        })
    }
    factor_feature <- colnames(data$fit$x)[sapply(data$fit$x, is.factor)]
    if(missing(feature))
        feature <- factor_feature
    if(any(!feature %in% factor_feature))
        stop(sprintf("%s are not factors in the data set.",
             example_string("feature", feature[!feature %in% factor_feature])))
    columns <- Map(convert, data$fit$x,
                   base[colnames(data$fit$x)], !colnames(data$fit$x) %in% names(base),
                   drop[colnames(data$fit$x)], !colnames(data$fit$x) %in% names(drop),
                   colnames(data$fit$x), MoreArgs=list(warn=TRUE))
    data$feature_selection <- rep(data$feature_selection,
        vapply(columns, function(x) if(!is.null(ncol(x))) ncol(x) else 1L, integer(1)))
    data$fit$x <- do.call(cbind, columns)
    names(data$feature_selection) <- colnames(data$fit$x)
    rm(columns)
    if(nrow(data$test$x) > 0){
        data$test$x <- do.call(cbind, 
            Map(convert, data$test$x,
                base[colnames(data$test$x)], !colnames(data$test$x) %in% names(base),
                drop[colnames(data$test$x)], !colnames(data$test$x) %in% names(drop),
                colnames(data$test$x), MoreArgs=list(warn=FALSE))
        )
    } else {
        data$test$x <- data$fit$x[FALSE,]
    }
    stopifnot(identical(colnames(data$fit$x), colnames(data$test$x)))
    data
}

#' Print log message during pre-processing
#' 
#' @param data Pre-processed data set.
#' @param ... Sent to \code{\link{log_message}}
#' @return The same data set as inputted. This only purpose of this function is
#'   to print a log message as a side effect.
#' @seealso \code{\link{pre_process}}, \code{\link{log_message}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pre_log_message <- function(data, ...){
    log_message(...)
    data
}

