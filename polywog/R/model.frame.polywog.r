##' Extracts the model frame from a fitted \code{\link{polywog}} model, as
##' \code{\link{model.frame.lm}} does for a fitted \code{\link{lm}} model.
##' @title Model frame of a polywog model
##' @param formula a fitted model of class \code{"polywog"} (the argument is
##' named \code{formula} for consistency with the generic function
##' \code{\link{model.frame}})
##' @param ... other arguments, currently ignored (but may later be adapted
##' for use as in \code{\link{model.frame.lm}})
##' @return A data frame containing the variables used to fit the model, with
##' additional attributes (e.g., \code{"terms"}) used to construct a model
##' matrix.
##' @seealso \code{\link{model.matrix.polywog}} for constructing the design
##' matrix.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method model.frame polywog
##' @export
model.frame.polywog <- function(formula, ...)
{
    ## Return the "model" element of the fitted model if it exists, otherwise
    ## try to reconstruct it from the environment in which the model was fit
    if (!is.null(formula$model)) {
        ans <- formula$model
    } else {
        ## Closely adapted from the code of 'model.frame.lm'
        fcall <- formula$call
        m <- match(c("data", "subset", "weights", "na.action"), names(fcall), 0L)
        fcall <- fcall[c(1L, m)]
        fcall$drop.unused.levels <- TRUE
        fcall[[1L]] <- quote(stats::model.frame)
        fcall$xlev <- formula$xlevels
        fcall$formula <- formula$formula
        env <- environment(formula$terms)
        if (is.null(env)) {
            env <- parent.frame()
        }
        ans <- eval(fcall, env, parent.frame())
    }

    ans
}
