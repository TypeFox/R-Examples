#' Setup a modeling procedure
#'
#' A modeling procedure is an object containing all information necessary to
#' carry out and evaluate the performance of a predictive modeling task with
#' \code{\link{fit}}, \code{\link{tune}}, or \code{\link{evaluate}}.
#' To use an out-of-the box algorithm with default values, only the
#' \code{method} argument needs to be set. See \code{\link{emil}} for a
#' list of available methods. To deviate from the defaults, e.g. by tuning
#' parameters or using a custom function for model fitting, set the appropriate
#' parameters as described below.
#' For a guide on how to implement a custom method see the documentaion page
#' \code{\link{extension}}.
#' 
#' @param method The name of the modeling method. Only needed to identify
#'   plug-in functions, i.e. if you supply them yourself there is no need to
#'   set \code{method}.
#' @param parameter A list of model parameters. These will be fed to the fitting
#'   function after the dataset (\code{x} and \code{y} parameters). To tune a
#'   parameter, supply the candidate values in a vector or list.
#' 
#'   When tuning more than one parameter, all combinations of parameter values
#'   will be tested, if the elements of \code{parameter} are named. To manually
#'   specify which parameter value combinations to try, leave the the elements
#'   unnamed (see example 3 and 4).
#'   
#'   Parameters that should have vectors or lists as values, e.g. \code{trControl}
#'   when using
#'   \code{\link{fit_caret}} to train pkg{caret} models, must be wrapped in an
#'   additional list. That is, to set a parameter value to a list, but not tune it,
#'   make it a list of length 1 containing the list to be used (see example 6).
#' @param fit_fun The function to be used for model fitting.
#' @param predict_fun The function to be used for model prediction.
#' @param importance_fun The function to be used for calculating or extracting 
#'   feature importances. See \code{\link{get_importance}} for details.
#' @param error_fun Performance measure used to evaluate procedures
#'   and to tune parameters. See \code{\link{error_fun}} for details.
#' @return An object of class \code{modeling_procedure}.
#' @example examples/modeling_procedure.r
#' @seealso \code{\link{emil}}, \code{\link{evaluate}},
#'   \code{\link{fit}}, \code{\link{tune}},
#'   \code{\link[=predict.model]{predict}}, \code{\link{get_importance}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
modeling_procedure <- function(method, parameter=list(), error_fun=NULL, fit_fun, predict_fun, importance_fun){
    if(any(sapply(list(NA, FALSE), identical, parameter))){
        warning("`parameter` must be supplied as a list. Assuming you really want `list()` i.e. to not set any parameters.")
        parameter <- list()
    }
    if(!is.null(names(parameter))){
        parameter <- apply(do.call(expand.grid, lapply(parameter, seq_along)), 1, function(i){
            Map("[[", parameter, i)
        })
    }
    if(missing(fit_fun)){
        fit_fun <- tryCatch({
            get(sprintf("fit_%s", method), globalenv())
        }, error = function(...){
            tryCatch({
                get(sprintf("fit_%s", method))
            }, error = function(err){
                err
            })
        })
    }
    if(!is.function(fit_fun))
        stop("No fitting function found.")
    if(missing(predict_fun)){
        predict_fun <- tryCatch({
            get(sprintf("predict_%s", method), globalenv())
        }, error = function(...){
            tryCatch({
                get(sprintf("predict_%s", method))
            }, error = function(err){
                err
            })
        })
    }
    if(missing(importance_fun)){
        importance_fun <- tryCatch({
            get(sprintf("importance_%s", method), globalenv())
        }, error = function(...){
            tryCatch({
                get(sprintf("importance_%s", method))
            }, error = function(err){
                err
            })
        })
    }
    structure(class = "modeling_procedure", .Data = list(
        method = if(missing(method)) "custom" else method,
        parameter = if(length(parameter) == 0) list() else if(length(parameter) == 1) parameter[[1]] else NULL,
        tuning = if(length(parameter) < 2) NULL else list(parameter = parameter, error = NULL),
        fit_fun = fit_fun,
        predict_fun = predict_fun,
        importance_fun = importance_fun,
        error_fun = error_fun
    ))
}

#' Print method for modeling procedure
#' 
#' @method print modeling_procedure
#' @param x modeling procedure.
#' @param ... Ignored (kept for S3 consistency).
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
print.modeling_procedure <- function(x, ...){
    yn <- function(x){
        if(is.null(x) || inherits(x, "error")){
            "no"
        } else {
            if(isdebugged(x)) "yes (debug)" else "yes"
        }
    }
    cat(sprintf(
"`%s` modeling procedure.

   model fitting function:       %s
   prediction function:          %s
   feature importance function:  %s
   individual error function:    %s
   
   number of parameter sets to tune over: %i
   tuned: %s\n",
        x$method,
        yn(x$fit_fun), yn(x$predict_fun), yn(x$importance_fun), yn(x$error_fun),
        if(is.null(x$tuning$parameter)) 1 else length(x$tuning$parameter),
        if(is_tunable(x)) if(is_tuned(x)) "yes" else "no" else "not needed"
    ))
}

#' Coerce to modeling procedure
#' 
#' @param x modeling procedure.
#' @param ... Ignored (kept for S3 consistency).
#' @return Modeling procedure
#' @examples
#' as.modeling_procedure("lda")
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
as.modeling_procedure <- function(x, ...){
    UseMethod("as.modeling_procedure")
}

#' @method as.modeling_procedure modeling_procedure
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
#' @noRd
as.modeling_procedure.modeling_procedure <- function(x, ...){
    x
}

#' @method as.modeling_procedure default
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
#' @noRd
as.modeling_procedure.default <- function(x, ...){
    stop("Cannot convert an object of class [", paste(class(x), collapse=","), "] to modeling_procedure.")
}

#' @method as.modeling_procedure character
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
#' @noRd
as.modeling_procedure.character <- function(x, ..., simplify=TRUE){
    if(length(x) == 1 && simplify){
        modeling_procedure(method = x)
    } else {
        if(is.null(names(x))) names(x) <- x
        lapply(x, modeling_procedure)
    }
}

#' Wrap a modeling procedures in a named list
#' 
#' The functions \code{\link{fit}}, \code{\link{tune}}, and
#' \code{\link{evaluate}} all work on lists of procedures rather than individual
#' procedures. This functions uniforms the structure of such lists from the
#' various form the user is allowed to input.
#' 
#' @param procedure Modeling procedure, as returned by
#'   \code{\link{modeling_procedure}}.
#' @return A named list of modeling procedures.
#' @noRd
#' @author Christofer \enc{Bäcklin}{Backlin}
multify <- function(procedure){
    multi_procedure <- (is.list(procedure) &&
                           !inherits(procedure, "modeling_procedure")) ||
                       (!is.list(procedure) && length(procedure) > 1)
    procedure <- if(inherits(procedure, "modeling_procedure")){
        list(procedure)
    } else {
        lapply(procedure, as.modeling_procedure)
    }
    name <- names(procedure)
    if(is.null(name)) name <- rep(NA, length(procedure))
    name <- ifelse(name %in% c(NA, ""), sapply(procedure, "[[", "method"), name)

    names(procedure) <- name
    attr(procedure, "multiple") <- multi_procedure
    procedure
}

#' Get debug flags of an object
#'
#' Normally you don't need this function when working with modeling procedures,
#' but in some special cases the flags are lost and need to be restored.
#'
#' @param x A complex object contianing functions.
#' @return A list of the same structure as \code{x} containing logicals
#'  indicating if each object is a function under debug.
#' @examples
#' m1 <- modeling_procedure("randomForest")
#' m2 <- modeling_procedure("pamr")
#' debug(m1$fit_fun)
#' isdebugged(m2$fit_fun)
#' debug_flag(m2) <- debug_flag(m1)
#' isdebugged(m2$fit_fun)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
debug_flag <- function(x){
    if(inherits(x, "modeling_procedure")){
        sapply(x, function(p) is.function(p) && isdebugged(p))
    } else if(is.list(x)){
        lapply(x, debug_flag)
    } else {
        stop("Invalid procedure.")
    }
}

#' @param value A list of logicals to be used to set debug flags.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
`debug_flag<-` <- function(x, value){
    if(inherits(x, "modeling_procedure") && is.logical(value)){
        for(f in names(which(value))) debug(x[[f]])
    } else if(is.list(x) && is.list(value)){
        for(i in seq_along(x)) debug_flag(x[[i]]) <- value[[i]]
    } else {
        stop("Object and debug flags doesn't match.")
    }
    x
}
