# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Model selection based on cross-validation
#' 
#' Combine cross-validation results for various models into one object and 
#' select the model with the best prediction performance.
#' 
#' Keep in mind that objects inheriting from class \code{"cv"} or 
#' \code{"cvSelect"} may contain multiple columns of cross-validation 
#' results.  This is the case if the response is univariate but the 
#' \code{\link[stats]{predict}} method of the fitted model returns a 
#' matrix.  
#' 
#' The \code{.reshape} argument determines how to handle such objects.  If 
#' \code{.reshape} is \code{FALSE}, all objects are required to have the same 
#' number of columns and the best model for each column is selected.  A typical 
#' use case for this behavior would be if the investigated models contain 
#' cross-validation results for a raw and a reweighted fit.  It might then be 
#' of interest to researchers to compare the best model for the raw estimators 
#' with the best model for the reweighted estimators.
#' 
#' If \code{.reshape} is \code{TRUE}, objects with more than one column of 
#' results are first transformed with \code{\link{cvReshape}} to have only one 
#' column.  Then the best overall model is selected.
#' 
#' It should also be noted that the argument names of \code{.reshape}, 
#' \code{.selectBest} and \code{.seFacor} start with a dot to avoid conflicts 
#' with the argument names used for the objects containing cross-validation 
#' results.
#' 
#' @param \dots  objects inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contain cross-validation results.
#' @param .reshape  a logical indicating whether objects with more than one 
#' column of cross-validation results should be reshaped to have only one 
#' column (see \dQuote{Details}).
#' @param .selectBest  a character string specifying a criterion for selecting 
#' the best model.  Possible values are \code{"min"} (the default) or 
#' \code{"hastie"}.  The former selects the model with the smallest prediction 
#' error.  The latter is useful for nested models or for models with a tuning 
#' parameter controlling the complexity of the model (e.g., penalized 
#' regression).  It selects the most parsimonious model whose prediction error 
#' is no larger than \code{.seFactor} standard errors above the prediction error 
#' of the best overall model.  Note that the models are thereby assumed to be 
#' ordered from the most parsimonious one to the most complex one.  In 
#' particular a one-standard-error rule is frequently applied.
#' @param .seFactor  a numeric value giving a multiplication factor of the 
#' standard error for the selection of the best model.  This is ignored if 
#' \code{.selectBest} is \code{"min"}.
#' 
#' @aliases print.cvSelect
#' 
#' @returnClass cvSelect
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer vector giving the number of folds used in 
#' cross-validation for the respective model.
#' @returnItem R  an integer vector giving the number of replications used in 
#' cross-validation for the respective model.
#' @returnItem best  an integer vector giving the indices of the models with 
#' the best prediction performance.
#' @returnItem cv  a data frame containing the estimated prediction errors for 
#' the models.  For models for which repeated cross-validation was performed, 
#' those are average values over all replications.
#' @returnItem se  a data frame containing the estimated standard errors of the 
#' prediction loss for the models.
#' @returnItem selectBest  a character string specifying the criterion used for 
#' selecting the best model.
#' @returnItem seFactor  a numeric value giving the multiplication factor of 
#' the standard error used for the selection of the best model.
#' @returnItem reps  a data frame containing the estimated prediction errors 
#' from all replications for those models for which repeated cross-validation 
#' was performed.  This is only returned if repeated cross-validation was 
#' performed for at least one of the models.
#' 
#' @note Even though the function allows to compare cross-validation results 
#' obtained with a different number of folds or a different number of 
#' replications, such comparisons should be made with care.  Hence warnings 
#' are issued in those cases.  For maximum comparability, the same data folds 
#' should be used in cross-validation for all models to be compared.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Hastie, T., Tibshirani, R. and Friedman, J. (2009) \emph{The Elements of 
#' Statistical Learning: Data Mining, Inference, and Prediction}.  Springer, 
#' 2nd edition.
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvTuning}}
#' 
#' @example inst/doc/examples/example-cvSelect.R
#' 
#' @keywords utilities
#' 
#' @export

cvSelect <- function(..., .reshape = FALSE, .selectBest = c("min", "hastie"), 
        .seFactor = 1) {
    ## initializations
    objects <- list(...)
    m <- length(objects)
    if(m == 0) stop("empty list of objects")
    # check class of objects
    isCvSelect <- sapply(objects, inherits, "cvSelect")
    if(!all(sapply(objects, inherits, "cv") | isCvSelect)) {
        stop("all objects must inherit from class \"cv\" or \"cvSelect\"")
    }
    .selectBest <- match.arg(.selectBest)
    # remove empty objects
    keep <- sapply(objects, function(x) ncv(x) > 0 && !isTRUE(nfits(x) == 0))
    objects <- objects[keep]
    m <- length(objects)
    if(m == 0) stop("all objects are empty")
    isCvSelect <- isCvSelect[keep]
    # check names for the supplied objects
    fits <- names(objects)
    if(is.null(fits)) {
        fits <- defaultFitNames(m)
    } else if(any(i <- fits == "")) fits[i] <- defaultFitNames(m)[i]
    names(objects) <- fits
    # check dimensions or reshape objects with more than one column
    d <- sapply(objects, ncv)
    if(isTRUE(.reshape)) {
        .reshape <- which(d > 1)
        objects[.reshape] <- lapply(objects[.reshape], cvReshape)
        isCvSelect[.reshape] <- TRUE
    } else if(length(unique(d)) > 1) {
        stop("all objects must have the same dimension")
    }
    ## prepare objects of class "cvSelect"
    if(any(isCvSelect)) {
        # prepare names
        fits <- as.list(fits)
        fits[isCvSelect] <- mapply(function(f, x) paste(f, x$cv$Fit, sep="."), 
            fits[isCvSelect], objects[isCvSelect], SIMPLIFY=FALSE)
        fits <- unlist(fits)
        # prepare basic information
        objects[isCvSelect] <- lapply(objects[isCvSelect], 
            function(x) {
                m <- nrow(x$cv)  # number of fits in current object
                x$n <- rep(x$n, length.out=m)
                x$K <- rep(x$K, length.out=m)
                x$R <- rep(x$R, length.out=m)
                # remove column specifying fit from results
                x$cv <- x$cv[, -1, drop=FALSE]
                x$se <- x$se[, -1, drop=FALSE]
                if(!is.null(x$reps)) x$reps <- x$reps[, -1, drop=FALSE]
                x
            })
    }
    ## combine basic information
    n <- unique(unlist(lapply(objects, function(x) x$n), use.names=FALSE))
    if(length(n) > 1) stop("different numbers of observations")
    K <- unlist(lapply(objects, function(x) x$K), use.names=FALSE)
    if(length(unique(K)) > 1) warning("different number of folds")
    R <- unlist(lapply(objects, function(x) x$R), use.names=FALSE)
    if(length(unique(R)) > 1) warning("different number of replications")
    names(K) <- names(R) <- fits
    ## combine CV results and standard errors
    cv <- lapply(objects, 
        function(x) {
            cv <- x$cv                                     # extract CV results
            if(is.null(dim(cv))) t(cv) else as.matrix(cv)  # return matrix
        })
    se <- lapply(objects, 
        function(x) {
            se <- x$se                                     # extract standard errors
            if(is.null(dim(se))) t(se) else as.matrix(se)  # return matrix
        })
    if(m > 1) {
        # check if names are the same for all objects
        cvNames <- colnames(cv[[1]])
        otherNames <- lapply(cv[-1], colnames)
        adjustNames <- !all(sapply(otherNames, identical, cvNames))
        if(adjustNames) cvNames <- defaultCvNames(length(cvNames))
    }
    cv <- do.call("rbind", cv)
    se <- do.call("rbind", se)
    if(m > 1 && adjustNames) {
        colnames(cv) <- colnames(se) <- cvNames
    }
    cv <- data.frame(Fit=factor(fits, levels=fits), cv, row.names=NULL)
    se <- data.frame(Fit=factor(fits, levels=fits), se, row.names=NULL)
    ## combine repeated CV results
    haveReps <- any(i <- sapply(objects, function(x) !is.null(x$reps)))
    if(haveReps) {
        # FIXME: safer solution that does not require the correct order of fits
        reps <- lapply(objects[i], function(x) as.matrix(x$reps))
        reps <- do.call("rbind", reps)
        if(m > 1 && adjustNames) colnames(reps) <- cvNames
        i <- which(R > 1)
        reps <- data.frame(Fit=factor(rep(fits[i], R[i]), levels=fits), 
            reps, row.names=NULL)
    }
    ## find best model
    if(.selectBest == "min") {
        .seFactor <- NA
        best <- sapply(cv[, -1, drop=FALSE], selectMin)
    } else {
        .seFactor <- rep(.seFactor, length.out=1)
        best <- sapply(names(cv)[-1], 
            function(j) selectHastie(cv[, j], se[, j], seFactor=.seFactor))
    }
    ## construct return object
    out <- list(n=n, K=K, R=R, best=best, cv=cv, se=se, 
        selectBest=.selectBest, seFactor=.seFactor)
    if(haveReps) out$reps <- reps
    class(out) <- "cvSelect"
    out
}
