##' Marginal effects for polywog models
##'
##' Computes average and observationwise marginal effects from a fitted
##' \code{\link{polywog}} model.
##'
##' For input variables that are binary, logical, or factors,
##' \code{margEff.polywog} computes a first difference with comparison to a
##' reference category.  All other variables are treated as continuous:
##' the function computes the partial derivative of the fitted value with
##' respect to the selected variable.
##' @param object a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.  The object must have a \code{model} element,
##' meaning it was fit with \code{model = TRUE}.
##' @param xvar a character string containing the name of a raw input variable
##' (from \code{object$varNames}).  Partial matches are allowed.
##' @param drop logical: whether to convert one-column matrices in the output to
##' vectors.
##' @param ... other arguments, currently ignored.
##' @return If \code{xvar} is specified, a numeric object containing
##' the marginal effect of the chosen variable at each observation in
##' \code{object$model}.  For factor variables, if there are more than two
##' levels or \code{drop = FALSE}, the returned object is a matrix; otherwise it
##' is a vector.
##'
##' If \code{xvar} is \code{NULL}, a list of such results for each raw input
##' variable in the model is returned.
##'
##' In either case, the returned object is of class \code{"margEff.polywog"}.
##' @seealso To plot the density of the observationwise marginal effects, see
##' \code{\link{plot.margEff.polywog}}.  For a table of average marginal effects
##' and order statistics, \code{\link{summary.margEff.polywog}}.
##'
##' To compute fitted values, see \code{\link{predict.polywog}} and
##' \code{\link{predVals}}.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method margEff polywog
##' @export
##' @example inst/examples/margEff.polywog.r
##' @importFrom miscTools margEff
margEff.polywog <- function(object, xvar = NULL, drop = FALSE, ...)
{
    ## If no variable specified, apply to each variable in the model
    if (is.null(xvar)) {
        ans <- lapply(object$varNames,
                      function(x) margEff(object, xvar = x, drop = drop))
        names(ans) <- object$varNames
        class(ans) <- "margEff.polywog"
        return(ans)
    }

    ## Extract original model data
    mf <- object$model
    if (is.null(mf))
        stop("Fitted object must contain 'model' to use margEff; re-run polywog with \"model = TRUE\"")

    ## Extract the variable of interest and determine its type
    xc <- getXcols(xvar, names(mf))
    xvar <- names(mf)[xc]  # replace with completed name
    x <- mf[, xc]
    xbinary <- is.numeric(x) && all(x %in% c(0, 1))

    if (is.numeric(x) && !xbinary) {
        ##-- CONTINUOUS VARIABLE --##

        X <- makeX(object$formula, mf)
        xc <- getXcols(xvar, colnames(X))
        ans <- computeMargEff(X = X,
                              poly_terms = object$polyTerms,
                              coef = coef(object),
                              coef_is_zero = coef(object) == 0,
                              xvar_col = xc - 1)

        ## Derivative of logit transformation
        if (object$family == "binomial") {
            xb <- predict(object, type = "link")
            ans <- dlogis(xb) * ans
        }
    } else {
        ##-- DISCRETE VARIABLE --##

        ## Switch to determine the "levels"
        if (is.factor(x)) {
            levs <- levels(x)
        } else if (is.character(x)) {
            warning("character variables should be converted to factors before fitting with polywog(); errors or unpredictable behavior may result")
            levs <- levels(as.factor(x))
        } else if (is.logical(x)) {
            levs <- c(TRUE, FALSE)
        } else if (xbinary) {
            levs <- c(1, 0)
        }

        ## Run predict() on the original data, setting the selected variable to
        ## each of its potential levels
        ans <- vector("list", length(levs))
        for (i in seq_along(levs)) {
            newmf <- mf
            newmf[, xc][] <- levs[i]

            ## Need to ensure that the "na.action" attribute of 'newdata' is
            ## NULL, or else padding may result if the model was originally run
            ## with na.exclude
            ans[[i]] <- predict(object,
                                newdata = structure(newmf, na.action = NULL),
                                type = "response")
        }
        ans <- do.call(cbind, ans)
        colnames(ans) <- levs

        ## Take first differences with respect to the reference category (last
        ## level)
        ans <- ans - ans[, ncol(ans)]
        ans <- ans[, -ncol(ans), drop = FALSE]  # Remove column of zeroes
        if (drop)
            ans <- base::drop(ans)
        attr(ans, "levels") <- levs
    }

    attr(ans, "xvar") <- xvar
    class(ans) <- "margEff.polywog"
    ans
}

## Convenience function to make a data frame out of a list of marginal effects
MEtoDF <- function(object)
{
    ## Figure out which entries are first differences
    whichFactors <- sapply(object, function(x) !is.null(attr(x, "levels")))

    ## Make intelligible names for first difference entries
    for (i in seq_len(sum(whichFactors))) {
        j <- which(whichFactors)[i]
        object[[j]] <- as.matrix(object[[j]])
        levs <- attr(object[[j]], "levels")
        xvar <- attr(object[[j]], "xvar")
        colnames(object[[j]]) <-
            paste(xvar, "=", levs[-length(levs)], " (vs. ",
                  levs[length(levs)], ")", sep = "")
    }

    as.data.frame(do.call(cbind, object))
}

##' Summarize marginal effects
##'
##' Generates a table of the average marginal effects and quartiles (or other
##' order statistics if requested) from a \code{"margEff.polywog"} object.
##' @param object output of \code{\link{margEff.polywog}}.
##' @param probs order statistics to display.
##' @param ... other arguments, currently ignored.
##' @return Table of results.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method summary margEff.polywog
##' @export
summary.margEff.polywog <- function(object, probs = seq(0, 1, by = 0.25), ...)
{
    if (!is.list(object)) {
        object <- list(object)
        names(object) <- attr(object[[1]], "xvar")
    }

    object <- MEtoDF(object)
    object <- sapply(object,
                     function(x) c("Mean" = mean(x), "SD" = sd(x),
                                   quantile(x, probs = probs)))
    object <- t(object)
    class(object) <- "summary.margEff.polywog"
    object
}

##' @S3method print summary.margEff.polywog
print.summary.margEff.polywog <- function(x,
                                          digits = max(3, getOption("digits") - 3),
                                          ...)
{
    printCoefmat(zapsmall(x), digits = digits, ...)
    invisible(x)
}
