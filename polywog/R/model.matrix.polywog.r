##' Constructs the design matrix used to fit a \code{\link{polywog}} model,
##' similar to \code{\link{model.matrix.lm}}.
##'
##' There are two types of model matrix a user might want to construct.
##' First, there is the matrix of the raw input terms that go into the
##' eventual polynomial expansion.  Such a matrix can be obtained by using
##' \code{type = "raw"} (the default). The other form of the model matrix is
##' the full polynomial expansion, where each column contains some power of
##' the raw inputs.  This can be obtained by using \code{type = "expanded"}.
##' @title Model matrix of a polywog model
##' @param object a fitted model of class \code{"polywog"}
##' @param type \code{"raw"}, the default, returns the non-expanded model
##' matrix with no intercept (same number of columns as
##' \code{object$polyTerms}).  \code{"expanded"} returns the polynomial
##' expansion used in fitting (number of columns equals
##' \code{length(object$coefficients)}).
##' @param ... other arguments to be passed to further methods (typically only
##' used internally)
##' @return The design matrix of the specified model, consisting either of raw
##' terms or the full polynomial expansion depending on the \code{type}
##' argument.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method model.matrix polywog
##' @export
model.matrix.polywog <- function(object, type = c("raw", "expanded"), ...)
{
    type <- match.arg(type)

    ## Use the stored X if available (i.e. polywog() was run with X=TRUE)
    ## Otherwise, reconstruct it from the model frame.  If polywog() was run
    ## with model=FALSE and the original data is no longer in memory,
    ## model.frame.polywog() will throw an error
    if (!is.null(object$X)) {
        X <- object$X
    } else {
        X <- makeX(object$formula, model.frame(object))
    }

    ## Compute the polynomial expansion if requested
    if (type == "expanded")
        X <- expandMatrix(X, object$polyTerms, intercept = TRUE)

    X
}
