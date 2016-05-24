#' PAMR adapted dataset pre-processing
#'
#' The predict framework is designed to work with dataset where rows correspond
#' to observations and columns to descriptors. PAMR wants it the other way, and
#' also to have the fitting set response vector supplied in a list with the
#' descriptors. This function applies a standard pre-processing function and
#' then reformats the result to satisfy PAMR.
#'
#' \code{pre_pamr} must be run last if chained with other pre-processing
#' functions, since it substantially reshapes the data.
#'
#' @param data Fitting and testing data sets, as returned by
#'   \code{\link{pre_split}}.
#' @return A list with fitting and testing sets, formatted the way pamr wants
#'   them.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{pre_process}}
#' @export
pre_pamr <- function(data){
    if(ncol(data$fit$x) == 1){
        notify_once(id = "pamr_univariate",
                    "PAMR not designed to handle univariate data. An all-zero dummy feature added as a workaround.",
                    fun = message)
    }
    data$fit$x <- structure(class = "pamr.data", .Data = list(
        x = if(ncol(data$fit$x) == 1){
            rbind(t(data$fit$x), dummy=0)
        } else {
            t(data$fit$x)
        },
        y = data$fit$y))
    data$test$x <- if(ncol(data$test$x) == 1){
        rbind(t(data$test$x), dummy=0)
    } else {
        t(data$test$x)
    }
    data
}

#' Fit nearest shrunken centroids model.
#'
#' Wrapped version of the \pkg{pamr} package implementation. Note that
#' this function uses internal cross-validation for determining the value
#' of the shrinkage threshold.
#'
#' @param x Dataset, numerical matrix with observations as rows.
#' @param y Class labels, factor.
#' @param error_fun Error function for tuning.
#' @param slim Set to \code{TRUE} if you want to return the fitted
#'   classifier but discard pamr's \code{cv.objects}, which can be large.
#'   memory efficient. This means that the element \code{cv$cv.objects} 
#'   containing the cross-validated fits will be dropped from the returned
#'   classifier.
#' @param cv Cross-validation scheme for shrinkage tuning. It should
#'   be supplied on one of the following forms:
#'   \itemize{
#'     \item{Resampling scheme produced with \code{\link{resample}}
#'       or \code{\link{resample_holdout}}.}
#'     \item{List with elements named \code{nrepeat} and \code{nfold}}
#'     \item{\code{NA}, \code{NULL} or \code{FALSE} to suppress shrinkage tuning.}
#'   }
#' @param nfold Sent to \code{\link[pamr]{pamr.cv}}. Only used if \code{cv} is missing.
#' @param threshold Shrinkage thresholds to try (referred to as 'lambda' in the
#'   literature). Chosen and tuned automatically by default, but must be given
#'   by the user if not tuned (see the \code{cv} argument) if you wish to use
#'   it with \code{\link{evaluate}}.
#' @param ... Sent to \code{\link[pamr]{pamr.train}}.
#' @param thres_fun Threshold selection function. Note that it is not uncommon
#'   that several thresholds will result in the same tuning error.
#' @return Fitted pamr classifier.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{predict_pamr}},
#'   \code{\link{importance_pamr}}, \code{\link{modeling_procedure}}
#' @export
fit_pamr <- function(x, y, error_fun, cv, nfold, threshold=NULL, ...,
                     thres_fun = function(thr, err) median(thr[err == min(err)]),
                     slim=FALSE){
    nice_require("pamr")
    if(!(is.list(x) && all(c("x", "y") %in% names(x)) && ncol(x$x) == length(x$y))){
        notify_once(id = "pamr_preprocess",
                    "Use the `pre_pamr` pre-processing function better pamr efficiency.",
                    fun = message)
        x <- list(x=t(x), y=y)
    }
    rm(y)
    if(missing(error_fun)){
        error_fun_frame <- sapply(sys.frames(), function(env){
            tryCatch({
                ef <- get("error_fun", envir=env)
                TRUE
            }, error=function(...) FALSE)
        })
        if(any(error_fun_frame)){
            error_fun <- get("error_fun", envir=sys.frames()[[max(which(error_fun_frame))]])
        } else {
            if(is.factor(x$y)){
                error_fun <- error_rate
            } else if(is.numeric(x$y)){
                error_fun <- rmse
            } else {
                stop("You must specify an error function!")
            }
        }
    }
    invisible(capture.output(
        tryCatch({
            model <- pamr::pamr.train(x, threshold=threshold, ...)
            if(missing(threshold) || length(threshold) > 1){
                if(missing(cv)){
                    model.cv <- if(missing(nfold)){
                        pamr::pamr.cv(model, x)
                    } else {
                        pamr::pamr.cv(model, x, nfold = nfold)
                    }
                } else if(is_blank(cv)){
                    stop("You cannot skip cross-validation when multiple thresholds are given.")
                } else {
                    if(!inherits(cv, c("crossvalidation", "holdout")))
                        cv <- resample("crossvalidation", x$y, nrepeat=cv$nrepeat, nfold=cv$nfold)
                    if(nrow(cv) != length(x$y))
                        stop("Resampling set for shrinkage selection does not match dataset in size.")
                    model.cv <- pamr::pamr.cv(model, x,
                        folds=lapply(cv, index_test))
                    if(slim){
                        model.cv$cv.objects <- NULL
                    }
                    model.cv$error <- sapply(seq_along(model.cv$threshold), function(i)
                        error_fun(model.cv$y, list(prediction=model.cv$yhat[[i]], probability=model.cv$probability[,,i])))
                }
            } else {
                if(!missing(cv) && !is_blank(cv))
                    notify_once(id = "pamr_ignoring_cv",
                                "Ignoring threshold tuning since only one threshold value was given.",
                                fun = message)
                model.cv <- NULL
            }
        }, error=function(...){
            stop(...)
        })
    ))
    return(list(model=model, cv=model.cv, thres_fun=thres_fun))
}


#' Prediction using nearest shrunken centroids.
#'
#' In case multiple thresholds give the same error the largest one is chosen
#' (i.e. the one keeping the fewest features).
#'
#' @param object Fitted classifier.
#' @param x Dataset of observations to be classified.
#' @param threshold Threshold to use for classification. This argument is only
#'   needed if you want to override the value set during model fitting.
#' @param thres_fun Threshold selection function. Only needed if you want to
#'   override the function set during model fitting.
#' @param ... Sent to \code{\link[pamr]{pamr.predict}}.
#' @return A list with elements:
#' \itemize{
#'     \item{\code{prediction}: Factor of predicted class memberships.}
#'     \item{\code{probability}: Data frame of predicted class probabilities.}
#' }
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_pamr}},
#'   \code{\link{importance_pamr}}, \code{\link{modeling_procedure}}
#' @export
predict_pamr <- function(object, x, threshold, thres_fun, ...){
    nice_require("pamr")
    if(nrow(x) != nrow(object$model$centroids)){
        if(ncol(x) != nrow(object$model$centroids))
            stop("PAMR takes datasets with observations as columns and descriptors as rows.")
        x <- t(x)
    }
    if(missing(threshold)){
        if(length(object$model$threshold) == 1){
            threshold <- object$model$threshold
        } else {
            if(missing(thres_fun)){
                thres_fun <- object$thres_fun
            }
            threshold <- thres_fun(object$cv$threshold, object$cv$error)
        }
    }
    list(prediction = pamr::pamr.predict(object$model, x, type="class", threshold=threshold, ...),
         probability = as.data.frame(pamr::pamr.predict(object$model,
                            x, type="posterior", threshold=threshold, ...)))
}


#' Feature importance of nearest shrunken centroids.
#' 
#' Calculated as the absolute difference between the overall centroid and a
#' class-wise shrunken centroid (which is the same for both classes except sign).
#'
#' In case multiple thresholds give the same error the largest one is chosen
#' (i.e. the one keeping the fewest features).
#' 
#' @param object Fitted pamr classifier
#' @param threshold Threshold to use for classification. This argument is only
#'   needed if you want to override the value set during model fitting.
#' @param thres_fun Threshold selection function. Only needed if you want to
#'   override the function set during model fitting.
#' @param ... Sent to \code{\link[pamr]{pamr.predict}}.
#' @return A data frame of feature importance scores.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_pamr}},
#'   \code{\link{predict_pamr}}, \code{\link{modeling_procedure}}
#' @export
importance_pamr <- function(object, threshold, thres_fun=max, ...){
    nice_require("pamr")
    if(missing(threshold)){
        if(length(object$model$threshold) == 1){
            threshold <- object$model$threshold
        } else {
            if(missing(thres_fun)){
                thres_fun <- object$thres_fun
            }
            threshold <- thres_fun(object$cv$threshold, object$cv$error)
        }
    }
    cen <- sweep(sweep(pamr::pamr.predict(object$model, , threshold, type="centroid", ...),
                       1, object$model$centroid.overall, "-"),
                 1, object$model$sd, "/")
    data.frame(feature = rownames(cen), cen, row.names=NULL)
}

