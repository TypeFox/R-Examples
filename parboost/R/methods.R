##' Print a parboost object
##'
##' Prints a basic description of a parboost object
##' @title Prints a short description of a parboost object.
##' @param x a parboost object.
##' @param ... Additional arguements passed to callies.
##' @return Prints a short descritpion of a parboost object.
##' @author Ronert Obst
##' @export
print.parboost <- function(x, ...) {

    cat("\n")
    cat("\t Distributed Boosting\n")
    cat("\n")
    if (!is.null(x$call))
        cat("Call to parboost:\n", deparse(x$call), "\n\n", sep = "")
    cat("Family: ", x$family, "\n")
    cat("Number of ensemble components", length(x$models), "\n")
    cat("Postprocessing method:", x$postprocessing, "\n")
    if (x$postprocessing == "none") {
        cat("Ensemble component weights:", x$ensemble_weights, "\n")
    } else {
        cat("Ensemble component weights (including intercept):", x$ensemble_weights, "\n")
    }
    cat("Split data:", x$split_data, "\n")
    cat("Number of folds used for cross validation:", x$folds, "\n")
    cat("Stepsize for tuning mstop:", x$stepsize_mstop, "\n")
    cat("Step size boosting: ", x$models[[1]]$control$nu, "\n")
    if (!is.null(x$models[[1]]$control$subsample)) {
        cat("Subsampled a fraction of", x$control$subsample, x$control$type, "in each iteration. \n")
    }
    cat("\n")
    invisible(x)
}

##' Print a summary of a parboost object
##'
##' Prints a basic summary of a parboost object
##' @title Prints a summary of a parboost object.
##' @param object a parboost object.
##' @param ... Additional arguements passed to callies.
##' @return Prints a summary of a parboost object.
##' @author Ronert Obst
##' @export
summary.parboost <- function(object, ...) {

    ret <- list(object = object, selprob = NULL)
    xs <- lapply(object$models, selected)
    nm <- lapply(object$models, variable.names)
    selprob <- lapply(1:length(object$models), function(i) tabulate(xs[[i]], nbins = length(nm[[i]])) / length(xs[[i]]))
    for (i in 1:length(object$models)) {
        names(selprob[[i]]) <- names(nm[[i]])
        selprob[[i]] <- sort(selprob[[i]], decreasing = TRUE)
    }

    ret$selprob <- lapply(1:length(object$models),
                          function(i) subset(selprob[[i]], selprob[[i]] > 0))
    class(ret) <- "summary.parboost"
    return(ret)
}

##' Print a summary of a parboost object
##'
##' Prints a basic summary of a parboost object
##' @title Prints a summary of a parboost object.
##' @param x a parboost object.
##' @param ... Additional arguements passed to callies.
##' @return Prints a summary of a parboost object.
##' @author Ronert Obst
##' @export
print.summary.parboost <- function(x, ...) {
    print(x$object)
    cat("Selection frequencies:\n")
    print(x$selprob)
    cat("\n")
}

##' Predict method for \code{parboost} objects
##'
##' If no new data is passed to \code{predict}, \code{predict} outputs
##' the fitted values. If you pass a data frame with new values to
##' \code{predict}, it will generate the predictions for them.
##' @title Generate predictions from parboost object
##' @param object Object of class parboost
##' @param newdata Optionally a data frame with new data to
##' predict
##' @param type String determining the type of prediction. The default
##' \code{"response"} is on the scale of the response variable and
##' \code{"link"} is on the scale of the predictors.
##' @param ... Additional parameters passed to predict.mboost.
##' @return Numeric vector of fitted values
##' @author Ronert Obst
##' @references T. Hothorn, P. Buehlmann, T. Kneib, M. Schmid, and
##' B. Hofner (2013). mboost: Model-Based Boosting, R package version
##' 2.2-3, \url{http://CRAN.R-project.org/package=mboost}.
##' @export
predict.parboost <- function(object, newdata=NULL, type = c("response", "link"), ...) {
    type <- match.arg(type)
    if(is.null(newdata)) {
        if (type == "link") {
            preds <- object$fitted
        } else {
            preds <- object$models[[1]]$family@response(object$fitted)
        }
        names(preds) <- NULL
        return(preds)
    } else {
        preds <- object$predict(object$models, newdata, type)
        names(preds) <- NULL
        return(preds)
    }
}

##' Extract coefficients from base learners that have a notion of coefficients
##'
##' Extract the coefficients of base learners which have a notion of
##' coefficients from each boosting model. Weighs the coefficients by
##' the postprocessed submodel weights.
##' @title Print coefficients for base learners with a notion of
##' coefficients
##' @param object Object of class parboost.
##' @param which Optionally a subset of base learners to evaluate as an
##' integer or character vector.
##' @param aggregate a character specifying how to aggregate
##' predictions or coefficients of single base learners. The default
##' returns the prediction or coefficient for the final number of
##' boosting iterations. \code{"cumsum"} returns a matrix with the
##' predictions for all iterations simultaneously (in
##' columns). \code{"none"} returns a list with matrices where the
##' \emph{j}th columns of the respective matrix contains the
##' predictions of the base learner of the \emph{j}th boosting
##' iteration (and zero if the base learner is not selected in this
##' iteration).
##' @param ... Additional arguements passed to callies.
##' @return Returns a list of coefficients
##' @author Ronert Obst
##' @references T. Hothorn, P. Buehlmann, T. Kneib, M. Schmid, and
##' B. Hofner (2013). mboost: Model-Based Boosting, R package version
##' 2.2-3, http://CRAN.R-project.org/package=mboost.
##' @method coef parboost
##' @export
coef.parboost <- function(object, which = NULL,
                          aggregate = c("sum", "cumsum", "none"), ...) {

    args <- list(...)
    if (length(args) > 0) {
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")
    }

    extract_coefs <- function(i, which, aggregate, weights, models) {
        if (object$postprocessing == "none") {
            coefs <- models[[i]]$coef(which = which, aggregate = aggregate)
            coefs <- lapply(coefs, `*`, weights[i])
        } else {
            coefs <- models[[i]]$coef(which = which, aggregate = aggregate)
            coefs <- lapply(coefs, `*`, weights[i+1])
        }
        return(coefs)
    }

    coefs_list <- lapply(1:length(object$models), extract_coefs, which, aggregate, object$ensemble_weights, object$models)
    coef_names <- unique(unlist(lapply(coefs_list, names)))

    coefs <- list()

    for (name in coef_names) {
        coefs[[name]] <- 0
        for (model in coefs_list) {
            for (cof in names(model)) {
                if (cof == name) {
                    coefs[[name]] <- coefs[[name]] + model[[name]]
                }
            }
        }
    }

    return(coefs)
}

##' Display selected base learners
##'
##' Displays the selected base learners from all submodels.
##' @title Selected base learners
##' @param object Object of class parboost
##' @param ... Parameters passed to selected.mboost
##' @return Numeric vector of selected base learners.
##' @author Ronert Obst
##' @references T. Hothorn, P. Buehlmann, T. Kneib, M. Schmid, and
##' B. Hofner (2013). mboost: Model-Based Boosting, R package version
##' 2.2-3, \url{http://CRAN.R-project.org/package=mboost}.
##' @export
selected.parboost <- function(object, ...) {
    unlist(lapply(object$models, selected))
}
