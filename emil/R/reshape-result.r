#' @importFrom tidyr extract_ gather_ spread_
NULL

#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Extractor function for modeling result
#' 
#' As opposed to the standard extractor function, this will keep the class.
#' 
#' @param x modeling result object, as produced by \code{\link{evaluate}}.
#' @param ... Sent to \link{Extract}.
#' @noRd
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
`[.modeling_result` <- function(x, ...){
    structure(unclass(x)[...], class="modeling_result")
}

#' Extract a subset of a tree of nested lists
#' 
#' Modeling results produced by \code{\link{evaluate}} comes in the
#' form of nested lists. This function can be used to subset or rearrange parts
#' of the results into vectors, matrices or data frames.
#' Also note the \code{\link[emil]{select}} function that provides an extension
#' to the \pkg{dplyr} package for data manipulation.
#' 
#' This function can only be used to extract data, not to assign.
#' 
#' @param x List of lists.
#' @param i Indexes to extract on the first level of the tree. Can also be a
#'   function that will be applied to the downstream result of the function.
#' @param ... Indexes to extract on subsequent levels.
#' @param error_value A template for the return value in case it is missing or
#'   invalid. Note that \code{NA} is a \code{\link{logical}} by default,
#'   causing \code{subtree} to also convert existing results to logicals.
#'   To get around this, please specify it as \code{as.numeric(NA)},
#'   \code{as.character(NA)}, or similar (see the example below).
#' @param warn Specifies whether warnings should be displayed (\code{0}),
#'   ignored (\code{-1}), or break execution (\code{1}). Works like the
#'   \code{\link{options}} parameter \code{warn}.
#' @param simplify Whether to collapse results into vectors or matrices when
#'   possible (\code{TRUE}) or to preserve the original tree structure as a
#'   list (\code{FALSE}).
#' @return A subset of the list tree.
#' @example examples/subtree.r
#' @seealso \code{\link{select}}, \code{\link{get_prediction}},
#'   \code{\link{get_importance}}, \code{\link{get_tuning}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
subtree <- function(x, i, ..., error_value, warn, simplify=TRUE){
    if(missing(error_value)) error_value <- NULL
    if(missing(warn)) warn <- is.null(error_value)
    if(is.null(error_value) && warn < 1){
        warning("With no `error_value` `warn` is ignored and all errors break the execution")
        warn <- 1
    }
    ret <- if(is.function(i)){
        if(missing(...)){
            i(x)
        } else {
            i(subtree(x, ..., error_value=error_value, warn=warn, simplify=simplify))
        }
    } else if(missing(...)){
        if(is.null(error_value)){
            x[i]
        } else {
            coerce_class <- function(x){
                x <- as(x, class(error_value))
                if(length(x) != length(error_value))
                    stop(sprintf("values must be length %i, but result is length %i",
                                 length(error_value), length(x)))
                x
            }
            lapply(x[i], function(xi){
                if(warn < 1) tryCatch({
                    coerce_class(xi)
                }, error = function(err){
                    if(warn == 0)
                        warning(err$message)
                    error_value
                }) else {
                    coerce_class(xi)
                }
            })
        }
    } else {
        lapply(x[i], subtree, ..., error_value=error_value, warn=warn, simplify=simplify)
    }
    if(simplify){
        if(length(ret) == 1){
            ret <- ret[[1]]
        } else if(all(sapply(ret, length) == 1)){
            ret <- unlist(ret, recursive=FALSE)
        } else if(all(sapply(lapply(ret, dim), is.null)) &&
                  all(sapply(ret, length) == length(ret[[1]]))){
            ret.class <- sapply(ret, class)
            ret.na <- sapply(ret, function(x) all(is.na(x)))
            i <- head(which(!ret.na), 1)
            if((is.numeric(ret[[i]]) || is.character(ret[[i]]) || is.logical(ret[[i]])) &&
                    length(unique(ret.class[!ret.na])) == 1){
                ret <- do.call(cbind, ret)
            }
        }
        if(is.null(dim(ret)) && !is.null(dim(x)) && length(ret) == length(x)){
            ret <- array(ret, dim=dim(x), dimnames=dimnames(x))
        }
    }
    ret
}

#' \pkg{emil} and \pkg{dplyr} integration
#' 
#' Modeling results can be converted to tabular format and manipulated using
#' \pkg{dplyr} and other Hadleyverse packages. This is accomplished by a class
#' specific \code{\link[dplyr]{select_}} function that differs somewhat in syntax
#' from the default \code{\link[dplyr]{select_}}.
#'
#' @param .data Modeling results, as returned by \code{\link{evaluate}}.
#' @param ... Not used, kept for consistency with \code{dplyr}.
#' @param .dots Indices to select on each level of \code{.data}, i.e.
#'   the first index specifies which top level elements of \code{.data} to
#'   select, the second specifies second-level-elements etc.
#'   The last index must select elements that can be converted to a data frame.
#'   In case the desired bottom-level element is related to the observations of
#'   a modeling task, e.g. the predctions of a test set, you must supply the
#'   resampling scheme used to produce \code{.data} at the appropriate level
#'   (see the examples).
#'
#'   The names of the \code{...} arguments specifies the names of the resulting
#'   data frame. Non-named arguments will be used to traverse the data but not
#'   returned.
#'
#'   In summary the \code{...} indices can be on the following forms:
#'   \describe{
#'     \item{Simple indices}{Anything that can be used to subset objects,
#'       e.g. integers, logicals, or characters.}
#'     \item{Functions}{A function that produces a data frame, vector or
#'       factor.}
#'     \item{Resampling schemes}{The same resampling scheme that was used to
#'       produce the modeling results.}
#'   }
#' @return A \code{\link{data.frame}} in long format.
#' @examples
#' # Produce some results
#' x <- iris[-5]
#' y <- iris$Species
#' names(y) <- sprintf("orchid%03i", seq_along(y))
#' cv <- resample("crossvalidation", y, nfold=3, nrepeat=2)
#' procedures <- list(nsc = modeling_procedure("pamr"),
#'                    rf = modeling_procedure("randomForest"))
#' result <- evaluate(procedures, x, y, resample=cv)
#' 
#' # Get the foldwise error for the NSC method
#' result %>% select(fold = TRUE, "nsc", error = "error")
#'
#' # Compare both methods 
#' require(tidyr)
#' result %>%
#'     select(fold = TRUE, method = TRUE, error = "error") %>%
#'     spread(method, error)
#' require(dplyr)
#' result %>%
#'     select(fold = TRUE, method = TRUE, error = "error") %>%
#'     group_by(method) %>% summarize(mean_error = mean(error))
#'
#' # Investigate the variability in estimated class 2 probability across folds
#' result %>%
#'     select(fold = cv, "nsc", "prediction", probability = function(x) x$probability[,2]) %>%
#'     spread(fold, probability)
#' @seealso subtree
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @name select
#' @importFrom dplyr select
#' @export
NULL

#' @rdname select
#' @method select_ list
#' @importFrom dplyr select_
#' @importFrom lazyeval lazy_eval
#' @export
select_.list <- function(.data, ..., .dots){
    select_list(.data, lazy_eval(.dots))
}
#' @method select_ modeling_result
#' @rdname select
#' @export
select_.modeling_result <- function(.data, ..., .dots){
    select_list(.data, lazy_eval(.dots))
}
#' @importFrom data.table rbindlist
select_list <- function(.data, .dots, id=NULL){
    named <- !is.null(names(.dots)) && names(.dots)[1] != ""
    subsetting <- inherits(.dots[[1]], c("numeric", "integer", "logical", "character"))

    if(length(.dots) == 1){
        if(subsetting){
            if(is.data.frame(.data) && (length(.dots[[1]]) > 1 || !named)){
                d <- .data[.dots[[1]]]
            } else {
                d <- data.frame(Value = .data[[.dots[[1]]]])
                if(named && length(d) > 0)
                    names(d) <- names(.dots)
            }
        } else if(is.function(.dots[[1]])){
            d <- .dots[[1]](.data)
            if(!inherits(d, c("data.frame", "NULL"))){
                if(named){
                    d <- data.frame(d)
                    names(d) <- names(.dots)
                } else {
                    stop("Functions must either be named or return a data frame or NULL.")
                }
            }
        } else {
            stop("Invalid argument.")
        }
        if(!is.null(id)){
            if(nrow(d) != length(id)) stop("Fold did not match data.")
            d$id <- id
        }
    } else {
        if(subsetting){
            d <- lapply(.data[.dots[[1]]], select_list, .dots[-1], id=id)
            if(named && length(d) > 0){
                if(is.null(names(d))){
                    d <- data.frame(..tmp = rep(seq_along(d), sapply(d, nrow)),
                                    rbindlist(d))
                } else {
                    d <- data.frame(..tmp = factor(rep(seq_along(d), sapply(d, nrow)),
                                              seq_along(d), names(d)),
                               rbindlist(d))
                }
                names(d)[1] <- names(.dots)[1]
            } else {
                d <- rbindlist(d)
            }
        } else if(is.function(.dots[[1]]) && length(.dots) > 1){
            stop("Not implemented.")
        } else if(inherits(.dots[[1]], "resample")){
            if(!identical(names(.data), names(.dots[[1]])))
                .data <- .data[names(.dots[[1]])]
            r <- if(identical(rownames(.dots[[1]]), as.character(1:nrow(.dots[[1]])))){
                seq_len(nrow(.dots[[1]]))
            } else {
                rownames(.dots[[1]])
            }
            d <- Map(function(x, i) select_list(x, .dots[-1], id=r[index_test(i)]),
                     .data, .dots[[1]])
            d <- data.frame(fold = factor(rep(names(.dots[[1]]), sapply(d, nrow)),
                                          names(.dots[[1]])),
                            rbindlist(d))
        } else {
            stop("Invalid argument.")
        }
    }
    d
}

#' Extract predictions from modeling results
#'
#' @param result Modeling result, as returned by \code{\link{evaluate}} and
#'   \code{\link{evaluate}}.
#' @param resample Resampling scheme used to create the results.
#' @param type The type of prediction to return. The possible types vary between
#'   modeling procedure.
#' @param format Table format of the output. See
#'   \url{http://en.wikipedia.org/wiki/Wide_and_narrow_data} for more info.
#' @return A data frame where the id column refers to the observations.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
get_prediction <- function(result, resample, type="prediction", format=c("long", "wide")){
    stopifnot(inherits(result, "modeling_result"))
    format <- match.arg(format)
    if(is_multi_procedure(result)){
        prediction <- select(result, fold=resample, method=TRUE, type,
                             prediction="prediction")
    } else {
        prediction <- select(result, fold=resample, type,
                             prediction=type)
        names(prediction)[names(prediction) == "prediction"] <- type
    }
    if(format == "wide"){
        spread_(prediction, "fold", type)
    } else {
        prediction
    }
}


#' Feature (variable) importance of a fitted model
#'
#' Note that different methods calculates feature importance in different
#' ways and that they are not directly comparable.
#'
#' When extending the \pkg{emil} framework with your own method, the importance
#' function should return a data frame where one column is called "feature" and
#' the remaining columns are named after the classes.
#'
#' @param object Fitted model.
#' @param format Table format of the output. See
#'   \url{http://en.wikipedia.org/wiki/Wide_and_narrow_data} for more info.
#' @param ... Sent on to the procedure's feature importance scoring function.
#' @return A vector of length p or an p-x-c matrix of feature importance
#'   scores where p is the number of descriptors and c is the number of classes.
#' @examples
#' procedure <- modeling_procedure("pamr")
#' model <- fit("pamr", x=iris[-5], y=iris$Species)
#' get_importance(model)
#' 
#' cv <- resample("crossvalidation", iris$Species, nrepeat=2, nfold=3)
#' result <- evaluate("pamr", iris[-5], iris$Species, resample=cv,
#'                    .save=c(importance=TRUE))
#' get_importance(result)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}
#' @export
get_importance <- function(object, format, ...){
    UseMethod("get_importance")
}
#' @method get_importance model
#' @export
get_importance.model <- function(object, format=c("wide", "long"), ...){
    format <- match.arg(format)
    imp <- object$procedure$importance_fun(object$model, ...)
    if(format == "long"){
        gather_(imp, "class", "importance", setdiff(colnames(imp), "feature"))
    } else {
        imp
    }
}
#' @method get_importance modeling_result
#' @export
get_importance.modeling_result <- function(object, format=c("wide", "long"), ...){
    format <- match.arg(format)
    reset_notification(id = "importance_missing")
    gatherer <- function(x){
        if(is.null(x$importance)){
            if(is.null(x$model)){
                stop("Could not extract importance since neither models nor importance was saved during evaluation.")
            }
            notify_once("importance_missing", 
                "Feature importance was not calculated during the evaluation. Calculating now.",
                fun=message)
            get_importance(x$model, "wide")
        } else {
            x$importance
        }
    }
    if(is_multi_procedure(object)){
        imp <- select(object, fold=TRUE, method=TRUE, gatherer)
    } else {
        imp <- select(object, fold=TRUE, gatherer)
    }
    if(format == "long"){
        gather_(imp, "class", "importance")
    } else {
        imp
    }
}
#' @method get_importance list
#' @export
get_importance.list <- function(object, format=c("wide", "long"), ...){
    format <- match.arg(format)
    lapply(object, get_importance, format=format, ...)
    # FIXME to one table
}

#' Extract parameter tuning statistics
#' 
#' @param object Fitted model or modeling procedure
#' @return A data frame of tuning statistics in long format.
#' @examples
#' procedure <- modeling_procedure("randomForest",
#'     parameter = list(mtry = c(1, 3),
#'                      nodesize = c(4, 10)))
#' model <- fit(procedure, x=iris[-5], y=iris$Species)
#' get_tuning(model)
#' 
#' options(emil_max_indent=4)
#' ho <- resample("holdout", iris$Species, nfold=5)
#' result <- evaluate(procedure, iris[-5], iris$Species, resample=ho,
#'                    .save=c(model=TRUE))
#' get_tuning(result)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
get_tuning <- function(object){
    UseMethod("get_tuning")
}
#' @method get_tuning model
#' @export
get_tuning.model <- function(object){
    get_tuning(object$procedure)
}
#' @method get_tuning modeling_procedure
#' @export
get_tuning.modeling_procedure <- function(object){
    stopifnot(is_tunable(object) && is_tuned(object))
    parameter_frame <- tryCatch({
        x <- Map(function(i, p) data.frame(parameter_set = i, as.data.frame(p)),
            seq_along(object$tuning$parameter),
            object$tuning$parameter)
        if(any(sapply(x, nrow) > 1))
            stop("Only scalar parameter values can be put in a tuning data frame.")
        do.call(rbind, x)
    }, error = function(...){
        warning("Could not convert parameter values to a data frame, using a numeric index instead.")
        data.frame(parameter_set = seq_along(object$tuning$parameter))
    })
    error_frame <- object$tuning$error %>%
            spread_("fold", "error")
    merge(parameter_frame, error_frame, by="parameter_set")
}
#' @method get_tuning modeling_result
#' @export
get_tuning.modeling_result <- function(object){
    tuning.folds <- names(object)
    if(is_multi_procedure(object)){
        method <- if(is.null(names(object[[1]]))){
            seq_along(object[[1]])
        } else {
            names(object[[1]])
        }
        names(method) <- method
        tuning <- lapply(method, function(m){
            parameter <- names(object[[1]][[m]]$model$procedure$tuning$parameter[[1]])
            tuning <- select(object, fold=TRUE, method=m, "model", get_tuning.model) %>%
                gather_("parameter", "value", parameter)
        })
        tuning <- do.call(rbind, tuning)
    } else {
        parameter <- names(object[[1]]$model$procedure$tuning$parameter[[1]])
        tuning <- select(object, fold=TRUE, "model", get_tuning.model) %>%
            gather_("parameter", "value", parameter)
    }
    tuning %>% gather_("tuning_fold", "error", tuning.folds)
}

#' Extract prediction performance
#' 
#' @param result Modeling result, as returned by \code{\link{evaluate}}.
#' @param format Table format of the output. See
#'   \url{http://en.wikipedia.org/wiki/Wide_and_narrow_data} for more info.
#' @return Data frame.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
get_performance <- function(result, format=c("wide", "long")){
    format <- match.arg(format)
    if(is_multi_procedure(result)){
        p <- select(result, fold=TRUE, method=TRUE, error="error")
        if(format == "wide"){
            spread_(p, "method", "error")
        } else {
            p
        }
    } else {
        select(result, fold=TRUE, error="error")
    }
}

