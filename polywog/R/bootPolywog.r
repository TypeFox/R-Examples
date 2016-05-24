##' Bootstrap a fitted polywog model
##'
##' Nonparametric bootstrap of the \code{\link{polywog}} regression procedure.
##' Can be run on a fitted model of class \code{"polywog"}, or within the
##' original procedure via the \code{boot} argument.  The function
##' \code{control.bp} can be used to pass options to \code{bootPolywog} when
##' bootstrapping within \code{\link{polywog}}.
##'
##' Parallel computation via the \code{.parallel} argument requires
##' registation of a backend for \code{\link[foreach:foreach]{\%dopar\%}}, as
##' in \code{\link{polywog}}.  In the case of \code{bootPolywog}, bootstrap
##' fitting is carried out in parallel, while cross-validation to choose the
##' penalization factor (assuming \code{reuse.lambda = FALSE}) is carried
##' out sequentially within each iteration.
##' @param model a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}} or the \code{"polywog.fit"} element of the output
##' of \code{\link{cv.polywog}}.
##' @param nboot number of bootstrap iterations.
##' @param .parallel logical: whether to perform computations in parallel
##' using a backend registered with \code{\link{foreach}}.
##' @param reuse.lambda logical: whether to use the penalization parameter from
##' the original fit (\code{TRUE}), or to cross-validate within each iteration
##' (\code{FALSE}, default).
##' @param reuse.penwt logical: whether to use the penalty weights from the
##' original fit (\code{TRUE}), or to re-calculate them within each iteration
##' (\code{FALSE}, default).
##' @param nlambda number of values of the penalty factor to examine in
##' cross-validation, as in \code{\link{polywog}}.
##' @param lambda.min.ratio ratio of the smallest value of the penalty factor
##' to the largest, as in \code{\link{polywog}}.
##' @param nfolds number of cross-validation folds to use.
##' @param thresh convergence threshold, as in \code{\link{polywog}}.  If
##' \code{NULL}, use the same value as in the original model.
##' @param maxit iteration limit for fitting, as in \code{\link{polywog}}.  If
##' \code{NULL}, use the same value as in the original model.
##' @param maxtries maximum number of attempts to generate a bootstrap sample
##' with a non-collinear model matrix (often problematic with lopsided binary
##' regressors) before stopping and issuing an error message.
##' @param min.prop for models with a binary response variable, minimum proportion of
##' non-modal outcome to ensure is included in each bootstrap iteration (for
##' example, set \code{min.prop = 0.1} to throw out any bootstrap iteration
##' where less than 10 percent or more than 90 percent of the observations are 1's).
##' @param report logical: whether to print a status bar.  Not available if
##' \code{.parallel = TRUE}.
##' @param .matrixOnly logical: whether to return just the matrix of bootstrap
##' coefficients (\code{TRUE}), or the originally supplied model with the
##' bootstrap matrix as the \code{boot.matrix} element (\code{FALSE}, default).
##' @return If \code{.matrixOnly = FALSE}, the returned object is \code{model}
##' with the bootstrap matrix included as its \code{boot.matrix} element.  If
##' \code{.matrixOnly = TRUE}, just the matrix is returned.  In either case, the
##' bootstrap matrix is a sparse matrix of class
##' \code{"\linkS4class{dgCMatrix}"}.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @example inst/examples/bootPolywog.r
##' @import foreach
##' @import iterators
##' @importFrom Matrix Matrix cBind
bootPolywog <- function(model,
                        nboot = 100,
                        .parallel = FALSE,
                        reuse.lambda = FALSE,
                        reuse.penwt = FALSE,
                        nlambda = 100,
                        lambda.min.ratio = 1e-4,
                        nfolds = 10,
                        thresh = NULL,
                        maxit = NULL,
                        maxtries = 1000,
                        min.prop = 0,
                        report = FALSE,
                        .matrixOnly = FALSE)
{
    ## Can't have a progress bar in parallel, sadly
    if (.parallel && report) {
        report <- FALSE
        warning("Cannot print status bar when parallelization is enabled")
    }

    ## Check for bad argument combinations
    if (model$method == "scad" && reuse.penwt)
        warning("Argument 'reuse.penwt' is ignored when method = \"scad\"")

    ## Obtain original model matrix and response variable
    X <- model.matrix(model, type = "raw")
    y <- model.response(model.frame(model))
    nobs <- model$nobs
    w <- model$weights
    if (is.null(w))
        w <- rep(1L, nobs)
    isBinary <- all(y %in% 0:1)

    ## Reuse lambda if told to, otherwise set it to NULL so that fitPolywog()
    ## will know to select it automatically via cross-validation
    lambda <- if (reuse.lambda) model$lambda else NULL

    ## If not specified, set convergence tolerance and maximum iterations to
    ## their defaults, depending on the type of model
    if (is.null(thresh))
        thresh <- model$thresh
    if (is.null(maxit))
        maxit <- model$maxit

    ## Set up progress bar
    pb <- txtProgressBar(min = 0, max = nboot)

    ## Use foreach() to iterate, potentially in parallel
    `%dofn%` <- if (.parallel) `%dopar%` else `%do%`
    ans <- foreach (i = icount(nboot)) %dofn% {
        tries <- 0

        ## Inner loop: keep drawing indices until we get a dataset that is
        ## linearly independent and passes the `min.prop` test (for binary
        ## outcomes), or until we hit `maxtries`
        repeat {
            tries <- tries + 1
            if (tries > maxtries) {
                stop("'maxtries' reached; no suitable bootstrap sample was found")
            }

            ind <- sample(seq_len(nobs), size = nobs, replace = TRUE)

            ## If the dependent variable is binary, check that it doesn't have
            ## too many 0s or 1s, as specified by the user
            meanY <- if (isBinary) mean(y[ind]) else NULL
            if (isBinary && (meanY < min.prop || 1 - meanY < min.prop)) {
                next
            }

            ## Compute the linear model coefficients and check for
            ## collinearity
            lmcoef <- tryCatch(computeLinearCoef(X = X[ind, , drop = FALSE],
                                                 y = y[ind],
                                                 polyTerms = model$polyTerms,
                                                 weights = w[ind],
                                                 allowRankDeficient = FALSE),
                               error = identity)
            if (!inherits(lmcoef, "error")) {
                break
            }
        }

        ## Reuse or compute penalty weights
        if (reuse.penwt) {
            penwt <- model$penwt
        } else {
            penwt <- computePenaltyWeights(X = X[ind, , drop = FALSE],
                                           y = y[ind],
                                           weights = w[ind],
                                           polyTerms = model$polyTerms,
                                           lmcoef = lmcoef,
                                           method = model$method,
                                           penwt.method = model$penwt.method,
                                           family = model$family,
                                           unpenalized = model$unpenalized)
        }

        ## Compute model fit
        fitPolywog <- switch(model$method,
                             alasso = fitAdaptiveLASSO,
                             scad = fitSCAD)
        fit <- fitPolywog(X = X[ind, , drop = FALSE],
                          y = y[ind],
                          weights = w[ind],
                          polyTerms = model$polyTerms,
                          family = model$family,
                          penwt = penwt,
                          lambda = lambda,
                          nlambda = nlambda,
                          lambda.min.ratio = lambda.min.ratio,
                          nfolds = nfolds,
                          foldid = NULL,
                          thresh = thresh,
                          maxit = maxit,
                          .parallel = FALSE)

        if (report)
            setTxtProgressBar(pb, i)

        Matrix(fit$coef, sparse = TRUE)
    }

    ## Print newline after progress bar completes
    if (report)
        cat("\n")

    ## Construct the bootstrap matrix, with each column being a set of
    ## coefficients (*not* each row as in v0.3.0 and earlier)
    ans <- do.call(cBind, ans)
    rownames(ans) <- names(coef(model))

    ## If .matrixOnly (typically only used within 'polywog'), return only the
    ## bootstrap matrix itself; otherwise, return the original model object with
    ## the bootstrap matrix included as an element
    if (.matrixOnly) {
        return(ans)
    } else {
        model$boot.matrix <- ans
        return(model)
    }
}

##' @rdname bootPolywog
##' @export
control.bp <- function(.parallel = FALSE,
                       reuse.lambda = FALSE,
                       reuse.penwt = FALSE,
                       nlambda = 100,
                       lambda.min.ratio = 1e-4,
                       nfolds = 10,
                       thresh = NULL,
                       maxit = NULL,
                       maxtries = 1000,
                       min.prop = 0,
                       report = FALSE)
{
    list(.parallel = .parallel,
         reuse.lambda = reuse.lambda,
         reuse.penwt = reuse.penwt,
         nlambda = nlambda,
         lambda.min.ratio = lambda.min.ratio,
         nfolds = nfolds,
         thresh = thresh,
         maxit = maxit,
         maxtries = maxtries,
         min.prop = min.prop,
         report = report)
}
