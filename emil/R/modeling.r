#' Evaluate a modeling procedure
#' 
#' This function performs the important task of evaluating the performance of
#' a modeling procedure with resampling, including tuning and pre-processing
#' to not bias the results by information leakage.
#'
#' @param procedure Modeling procedure, or list of modeling procedures, as
#'   produced by \code{\link{modeling_procedure}}.
#' @param x Dataset, observations as rows and descriptors as columns.
#' @param y Response vector.
#' @param resample The test subsets used for parameter tuning. Leave blank to
#'   randomly generate a resampling scheme of the same kind as is used by
#'   \code{\link{evaluate}} to assess the performance of the whole
#'   modeling_procedure.
#' @param pre_process Function that performs pre-processing and splits dataset
#'   into fitting and test subsets.
#' @param .save What parts of the modeling results to return to the user. If
#'   \code{importance} is \code{FALSE} varible importance calculation will be
#'   skipped.
#' @param .cores Number of CPU-cores to use for parallel computation.
#'   The current implementation is based on \code{\link{mcMap}}, which
#'   unfortunatelly do not work on Windows systems. It can however be
#'   re-implemented by the user fairly easily by setting up a PSOCK cluster and
#'   calling \code{\link{parLapply}} as in the example below. This solution
#'   might be included in future versions of the package, after further
#'   investigation.
#' @param .checkpoint_dir Directory to save intermediate results to, after
#'   every completed fold. The directory will be created if it doesn't exist,
#'   but not recursively.
#' @param .return_error If \code{FALSE} the entire modeling is aborted upon an
#'   error. If \code{TRUE} the modeling of the particular fold is aborted and
#'   the error message is returned instead of its results.
#' @param .verbose Whether to print an activity log.
#' @return A list tree where the top level corresponds to folds (in case of
#'   multiple folds), the next level corresponds to the modeling procedures
#'   (in case of multiple procedures), and the final level is specified by the
#'   \code{.save} parameter. It typically contains a subset of the following
#'   elements:
#'   \describe{
#'       \item{\code{error}}{Performance estimate of the fitted model. See
#'           \code{\link{error_fun}} for more information.}
#'       \item{\code{fit}}{Fitted model.}
#'       \item{\code{prediction}}{Predictions given by the model.}
#'       \item{\code{importance}}{Feature importance scores.}
#'       \item{\code{tune}}{Results from the parameter tuning. See
#'           \code{\link{tune}} for details.}
#'   }
#' @example examples/evaluate.r
#' @seealso \code{\link{emil}}, \code{\link{modeling_procedure}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
evaluate <- function(procedure, x, y, resample, pre_process=pre_split,
    .save=c(model=TRUE, prediction=TRUE, error=TRUE, importance=FALSE),
    .cores=1, .checkpoint_dir=NULL, .return_error=.cores > 1,
    .verbose=getOption("emil_verbose", TRUE)){

    log_message(.verbose, "Evaluating modeling performance...")

    procedure <- multify(procedure)
    
    # Get the real response vector since we will need it to prepare the analysis
    y_original <- y
    y <- get_response(x, y)

    if(missing(resample))
        resample <- emil::resample("crossvalidation", y, nfold=2, nrepeat=2)
    if(inherits(resample, "fold")){
        multi.fold <- FALSE
        resample <- data.frame(resample)
        names(resample) <- attr(resample[[1]], "fold.name")
    } else {
        multi.fold <- TRUE
    }
    make.na <- is.na(y) & !Reduce("&", lapply(resample, is.na))
    if(any(make.na)){
        log_message(indent(.verbose, 1), "%i observations will be excluded from the modeling due to missing values.", sum(make.na))
        resample[make.na,] <- NA
    }

    .save.default <- c(model=FALSE, prediction=TRUE, error=TRUE, importance=FALSE)
    for(s in names(.save)){
        .save.default[pmatch(s, names(.save.default))] <- .save[s]
    }
    .save <- .save.default
    rm(.save.default)

#------------------------------------------------------------------------------o
#   Make sure all plug-ins exist before we start crunching

    do.tuning <- !sapply(procedure, is_tuned)
    missing.fun <- unlist(Map(function(p, to.be.tuned) c(
        if(!is.function(p$fit_fun))
            sprintf("fit_%s", p$method)
        else NULL,
        if(!is.function(p$predict_fun) && (.save["prediction"] || to.be.tuned))
            sprintf("predict_%s", p$method)
        else NULL,
        if(!is.function(p$importance_fun) && .save["importance"])
            sprintf("importance_%s", p$method)
        else NULL
    ), procedure, do.tuning))
    if(!is.null(missing.fun))
        stop(sprintf("Plug-in%s function %s not found.",
            if(length(missing.fun) > 1) "s" else "",
            paste("`", missing.fun, "`", sep="", collapse=", ")))

    for(i in seq_along(procedure)){
        if(is.null(procedure[[i]]$error_fun)){
            procedure[[i]]$error_fun <- if(class(y) == "factor"){
                error_rate
            } else if(inherits(y, "Surv")){
                neg_harrell_c
            } else if(is.numeric(y)){
                rmse
            } else {
                stop("Unknown type of response vector, cannot guess performance measure. Please set `error_fun` manually.")
            }
        }
    }

#------------------------------------------------------------------------------o
#   Set up parallelization, error handling and checkpointing

    if(.Platform$OS.type == "windows" && .cores > 1){
        notify_once(id = "windows_parallelization",
                    "Parallelization is not yet implemented for windows systems. Please set it up manually as described in `?evaluate`. Running sequentially.",
                    fun = warning)
        .cores <- 1
    }
    if(.cores > 1){
        requireNamespace("parallel")
        Map_fun <- function(f, ...)
            parallel::mcmapply(f, ..., SIMPLIFY=FALSE,
                               mc.silent=.verbose, mc.cores=.cores)
    } else {
        Map_fun <- Map
    }
    if(.return_error){
        modeling_fun <- function(f, ...){
            Map_fun(function(...){
                tryCatch(f(...), error = function(err){
                    log_message(indent(.verbose, 2), "An error occurred, skipping fold.")
                    err
                })
            }, ...)
        }
    } else {
        modeling_fun <- Map_fun
    }
    if(!is.null(.checkpoint_dir)){
        if(!file.exists(.checkpoint_dir))
            dir.create(.checkpoint_dir)
        if(file.access(.checkpoint_dir, 6) != 0)
            stop("You do not have read and write permissions to the checkpoint directory.")
        checkpoint.files <- sprintf("%s/%s.Rdata",
            .checkpoint_dir, gsub("\\W+", "-", names(resample)))
    } else {
        checkpoint.files <- list(NULL)
    }

#------------------------------------------------------------------------------o
#   Build and test models

    counter <- 0
    res <- structure(class="modeling_result", .Data=modeling_fun(function(fold, fold.name, checkpoint.file){

        # Setup run time estimation
        if(.cores == 1 && is.null(checkpoint.file)){
            counter <<- counter + 1
            if(counter == 1) t1 <- Sys.time()
        }

        # Print status message
        fold.message <- if(inherits(fold, "crossvalidation")){
            sub("^rep(\\d+)fold(\\d+)$", "Repeat \\1, fold \\2", fold.name)
        } else if(inherits(resample, "holdout")){
            sub("^fold(\\d+)$", "Fold \\1", fold.name)
        } else {
            fold.name
        }

        # Check for checkpoint files
        if(!is.null(checkpoint.file) && file.exists(checkpoint.file)){
            fold.message <- paste(fold.message, "already completed.")
            en <- new.env()
            load(checkpoint.file, envir=en)
            return(en$res)
        }
        log_message(indent(.verbose, 1), fold.message)

        # Disable further messages if run in parallel
        if(.cores > 1) .verbose <- 0

        # Do the work
        if(any(do.tuning)){
            # Tune before extracting the training and testing sets (below) to
            # not produce any unnecessary in-memory copies of the dataset.
            tune.subset <- subresample(fold = fold, y = y)
            fold.procedure <- tune(procedure, x, y, resample=tune.subset,
                pre_process=pre_process, .save=NULL, .verbose=indent(.verbose, 2))
        } else {
            fold.procedure <- procedure
        }
        log_message(indent(.verbose, 2), "Extracting fitting and testing datasets.")
        if(is.function(pre_process)){
            if(".verbose" %in% names(formals(pre_process))){
                sets <- pre_process(x, y_original, fold, .verbose = indent(.verbose, 2))
            } else {
                sets <- pre_process(x, y_original, fold)
            }
        } else if(is.list(pre_process)){
            sets <- pre_process[[1]](x, y, fold)
            if(length(pre_process) > 1) for(i in 2:length(pre_process)){
                sets <- pre_process[[i]](sets)
            }
            sets
        } else {
            stop("Invalid pre-processing.")
        }

        res <- lapply(seq_along(fold.procedure), function(i){
            # Rather than doing all models at once we do one at a time in case
            # they require a lot of memory. The [[-workaround is to keep the
            # procedure names printing in the correct way.
            model <- fit(fold.procedure[i], x=sets$fit$x, y=sets$fit$y,
                         .verbose=indent(.verbose, 2))[[1]]
            prediction <- predict(object=model, x=sets$test$x, .verbose=indent(.verbose, 1))
            list(error = if(.save["error"]){
                    fold.procedure[[i]]$error_fun(sets$test$y, prediction)
                 } else NULL,
                 model = if(.save["model"]) model else NULL,
                 prediction = if(.save["prediction"]) prediction else NULL,
                 importance = if(.save["importance"]) get_importance(model)
             )
        })
        names(res) <- names(fold.procedure)

        # Estimate run time
        if(.cores == 1 && is.null(checkpoint.file)){
            if(.verbose && counter == 1 && is.data.frame(resample) && ncol(resample) > 1){
                # Fold result size
                os <- object.size(res)
                os.i <- trunc(log(os)/log(1024))
                if(os.i == 0){
                    log_message(indent(.verbose, 1), "Result size is %i B.", os, time=FALSE)
                } else {
                    log_message(indent(.verbose, 1), "Result size is %.2f %s.",
                        exp(log(os) - os.i * log(1024)), c("KiB", "MiB", "GiB", "TiB", "PiB", "EiB")[os.i],
                        time=FALSE)
                }
                # Estimated completion time
                t2 <- t1 + difftime(Sys.time(), t1, units="sec")*ncol(resample)
                fmt <- if(difftime(t2, t1, units="days") < 1){
                    "%H:%M"
                } else if(difftime(t2, t1, units="days") < 2){
                    "%H:%M tomorrow"
                } else if(difftime(t2, t1, units="days") < 365){
                    "%H:%M, %b %d"
                } else {
                    "%H:%M, %b %d, %Y"
                }
                log_message(indent(.verbose, 1),
                          "Estimated completion time %sis %s.",
                          if("tune" %in% sapply(sys.calls(), function(x) as.character(x)[1]))
                              "of tuning " else "",
                          strftime(t2, fmt), time=FALSE)
            }
        }

        # Save checkpoint file
        if(!attr(procedure, "multiple")) res <- res[[1]]
        if(!is.null(checkpoint.file)){
            save(res, file=checkpoint.file)
        }

        res
    }, resample, names(resample), checkpoint.files))
    if(multi.fold) res else res[[1]]
}


#' Fit a model
#'
#' Fits a model according to a modeling procedure. If the procedure contains
#' untuned parameters they will automatically be tuned prior to fitting.
#' 
#' @param procedure Modeling procedure, or list of modeling procedures, as
#'   produced by \code{\link{modeling_procedure}}.
#' @param x Dataset, observations as rows and descriptors as columns.
#' @param y Response vector.
#' @param ... Sent to \code{\link{tune}}, in case tuning is required,
#'   which will pass them on to \code{\link{evaluate}}.
#' @param .verbose Whether to print an activity log. Set to \code{-1} to
#'   suppress all messages.
#' @return A list of fitted models.
#' @examples
#' mod <- fit("lda", x=iris[-5], y=iris$Species)
#' @seealso \code{\link{emil}}, \code{\link{modeling_procedure}},
#'   \code{\link{evaluate}}, \code{\link{tune}},
#'   \code{\link[=predict.model]{predict}}, \code{\link{get_importance}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
fit <- function(procedure, x, y, ..., .verbose=getOption("emil_verbose", FALSE)){
    #log_message(.verbose, "Fitting models...")
    procedure <- multify(procedure)

    missing.fun <- unlist(lapply(procedure, function(p)
        if(!is.function(p$fit_fun)) sprintf("fit_%s", p$method) else NULL))
    if(!is.null(missing.fun))
        stop(sprintf("Plug-in%s function %s not found.",
            if(length(missing.fun) > 1) "s" else "",
            paste("`", missing.fun, "`", sep="", collapse=", ")))

    need.tuning <- !sapply(procedure, is_tuned)
    if(any(need.tuning)){
        procedure[need.tuning] <- tune(procedure[need.tuning], x=x, y=y, ...,
            .verbose=indent(.verbose, 1))
    }
    res <- Map(function(p, mn){
        log_message(.verbose, "Fitting %s", mn)
        if(".verbose" %in% names(formals(p$fit_fun)))
            p$parameter$.verbose <- indent(.verbose, 1)
        model <- do.call(function(...) p$fit_fun(x=x, y=y, ...), p$parameter)
        structure(list(model = model, procedure = p), class="model")
    }, procedure, names(procedure))
    # Restore debug flags (for some reason Map removes them)
    for(i in seq_along(procedure))
        debug_flag(res[[i]]$procedure) <- debug_flag(procedure[[i]])
    if(attr(procedure, "multiple")) res else res[[1]]
}


#' Tune parameters of modeling procedures
#'
#' These functions are rarely needed to be called manually as they are
#' automatically called by \code{\link{fit}} and \code{\link{evaluate}}
#' when needed.
#' 
#' @param procedure Modeling procedure, or list of modeling procedures, as
#'   produced by \code{\link{modeling_procedure}}.
#' @param ... Sent to \code{\link{evaluate}}.
#' @param .verbose Whether to print an activity log. Set to \code{-1} to
#'   suppress all messages.
#' @return A tuned modeling procedures or a list of such.
#' @examples
#' procedure <- modeling_procedure("randomForest", parameter=list(mtry=1:4))
#' tuned.procedure <- tune(procedure, x=iris[-5], y=iris$Species)
#' mod <- fit(tuned.procedure, x=iris[-5], y=iris$Species)
#' @seealso \code{\link{emil}}, \code{\link{modeling_procedure}},
#'   \code{\link{evaluate}}, \code{\link{fit}},
#'   \code{\link[=predict.model]{predict}}, \code{\link{get_importance}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @importFrom dplyr group_by_ filter_ mutate_ sample_n
#' @export
tune <- function(procedure, ..., .verbose=getOption("emil_verbose", FALSE)){
    log_message(.verbose, "Tuning parameters...")
    procedure <- multify(procedure)

    do.tuning <- sapply(procedure, is_tunable)
    discard.tuning <- do.tuning & sapply(procedure, is_tuned)
    for(i in which(discard.tuning))
        procedure[[i]]$tuning <- NULL
    tune.procedure <- lapply(procedure[do.tuning], function(p){
        lapply(p$tuning$parameter, function(pp){
            p$parameter <- pp
            p["tuning"] <- list(NULL)
            p
        })
    })
    tune.name <- rep(names(tune.procedure), sapply(tune.procedure, length))
    tune.procedure <- unlist(tune.procedure, recursive=FALSE)
    names(tune.procedure) <- ave(tune.name, tune.name, FUN=function(x){
        if(length(x) == 1) x else sprintf("%s (%i)", x, seq_along(x))
    })

    tuning <- evaluate(tune.procedure, ..., .verbose=indent(.verbose, 1))
    procedure.id <- rep(which(do.tuning),
                   sapply(procedure[do.tuning], function(p) length(p$tuning$parameter)))
    for(i in which(do.tuning)){
        procedure[[i]]$tuning$error <-
            select(tuning, fold=TRUE, method = procedure.id == i, error = "error") %>%
            mutate_(parameter_set = "as.integer(method)") %>%
            select_("-method")
        best_parameter <- procedure[[i]]$tuning$error %>%
            group_by_("parameter_set") %>%
            summarize_(mean_error = "mean(error)") %>%
            filter_("mean_error == min(mean_error)") %>%
            sample_n(1)
        procedure[[i]]$parameter <-
            procedure[[i]]$tuning$parameter[[best_parameter$parameter_set]]
        procedure[[i]]$tuning$result <- lapply(tuning, "[", procedure.id == i)
    }
    if(attr(procedure, "multiple")) procedure else procedure[[1]]
}
#' @rdname tune
#' @return Logical indicating if the procedure(s) are tuned.
#' @export
is_tuned <- function(procedure){
    stopifnot(inherits(procedure, "modeling_procedure"))
    !is_tunable(procedure) || !is.null(procedure$parameter)
}
#' @rdname tune
#' @return Logical indicating if the has tunable parameters.
#' @export
is_tunable <- function(procedure){
    stopifnot(inherits(procedure, "modeling_procedure"))
    !is.null(procedure$tuning)
}
#' @rdname tune
#' @return A list of untuned modeling procedures.
#' @export
detune <- function(procedure){
    stopifnot(inherits(procedure, "modeling_procedure"))
    debug.flags <- sapply(procedure, function(p) is.function(p) && isdebugged(p))
    procedure$tuning <- NULL
    lapply(procedure[debug.flags], debug)
    procedure
}


#' Predict the response of unknown observations
#'
#' @method predict model
#' @param object Fitted model.
#' @param x Data set with observations whose response is to be predicted.
#' @param ... Sent to the procedure's prediction function.
#' @param .verbose Whether to print an activity log.
#' @return See the documentation of procedure's method.
#' @examples
#' mod <- fit("lda", x=iris[-5], y=iris$Species)
#' prediction <- predict(mod, iris[-5])
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{modeling_procedure}},
#'   \code{\link{evaluate}},\code{\link{fit}}, \code{\link{tune}},
#'   \code{\link{get_importance}}
#' @export
predict.model <- function(object, x, ..., .verbose=FALSE){
    if(".verbose" %in% names(formals(object$procedure$predict_fun))){
        p <- object$procedure$predict_fun(object = object$model, x = x, ..., .verbose=.verbose)
    } else {
        p <- object$procedure$predict_fun(object = object$model, x = x, ...)
    }
    structure(p, class="prediction")
}


#' Extract the response from a data set
#' 
#' @param x Data set features.
#' @param y Response vector or any other type of objects that describe how to
#'   extract the response vector from \code{x}.
#' @return A response vector.
#' @examples
#' identical(iris$Species, get_response(iris, "Species"))
#' identical(iris$Sepal.Length, get_response(iris, Sepal.Length ~ .))
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
get_response <- function(x, y){
    UseMethod("get_response", y)
}
#' @method get_response default
#' @export
get_response.default <- function(x, y){
    y
}
#' @method get_response default
#' @export
get_response.character <- function(x, y){
    x[,y]
}
#' @method get_response default
#' @export
get_response.formula <- function(x, y){
    model.response(model.frame(y, x))
}

