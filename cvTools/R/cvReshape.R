# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Reshape cross-validation results
#' 
#' Reshape cross-validation results into an object of class \code{"cvSelect"}
#' with only one column of results.
#' 
#' @aliases cvReshape.cv cvReshape.cvSelect
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results.
#' @param selectBest  a character string specifying a criterion for selecting 
#' the best model.  Possible values are \code{"min"} (the default) or 
#' \code{"hastie"}.  The former selects the model with the smallest prediction 
#' error.  The latter is useful for nested models or for models with a tuning 
#' parameter controlling the complexity of the model (e.g., penalized 
#' regression).  It selects the most parsimonious model whose prediction error 
#' is no larger than \code{seFactor} standard errors above the prediction error 
#' of the best overall model.  Note that the models are thereby assumed to be 
#' ordered from the most parsimonious one to the most complex one.  In 
#' particular a one-standard-error rule is frequently applied.
#' @param seFactor  a numeric value giving a multiplication factor of the 
#' standard error for the selection of the best model.  This is ignored if 
#' \code{selectBest} is \code{"min"}.
#' @param \dots  additional arguments to be passed down.
#' 
#' @returnClass cvSelect
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds used in 
#' cross-validation.
#' @returnItem R  an integer giving the number of replications used in 
#' cross-validation.
#' @returnItem best  an integer giving the index of the model with the best 
#' prediction performance.
#' @returnItem cv  a data frame containing the estimated prediction errors for 
#' the models.  For repeated cross-validation, those are average values over 
#' all replications.
#' @returnItem se  a data frame containing the estimated standard errors of the 
#' prediction loss for the models.
#' @returnItem selectBest  a character string specifying the criterion used for 
#' selecting the best model.
#' @returnItem seFactor  a numeric value giving the multiplication factor of 
#' the standard error used for the selection of the best model.
#' @returnItem reps  a data frame containing the estimated prediction errors 
#' for the models from all replications.  This is only returned if repeated 
#' cross-validation was performed.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Hastie, T., Tibshirani, R. and Friedman, J. (2009) \emph{The Elements of 
#' Statistical Learning: Data Mining, Inference, and Prediction}.  Springer, 
#' 2nd edition.
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, \code{\link{cvTuning}}
#' 
#' @example inst/doc/examples/example-cvReshape.R
#' 
#' @keywords utilities
#' 
#' @export

cvReshape <- function(x, ...) UseMethod("cvReshape")

#' @rdname cvReshape
#' @method cvReshape cv
#' @export
cvReshape.cv <- function(x, 
        selectBest = c("min", "hastie"), seFactor = 1, ...) {
    # initializations
    if(ncv(x) == 0 || isTRUE(nfits(x) == 0)) stop("empty object")
    cvNames <- cvNames(x)
    # create list of objects with one column
    cvName <- defaultCvNames(1)
    objects <- lapply(cvNames, 
        function(s) {
            xs <- subset(x, select=s)
            cvNames(xs) <- cvName
            xs
        })
    # substitute "CV" in default names by "Fit"
    if(identical(cvNames, defaultCvNames(length(cvNames)))) {
        fitName <- defaultFitNames(1)
        cvNames <- gsub(cvName, fitName, cvNames, fixed=TRUE)
    }
    # call cvSelect() to combine the model fits
    names(objects) <- cvNames
    objects$.selectBest <- selectBest
    objects$.seFactor <- seFactor
    out <- do.call(cvSelect, objects)
    if(length(out$K) > 1) out$K <- out$K[1]
    if(length(out$R) > 1) out$R <- out$R[1]
    out
}

#' @rdname cvReshape
#' @method cvReshape cvSelect
#' @export
cvReshape.cvSelect <- cvReshape.cv
