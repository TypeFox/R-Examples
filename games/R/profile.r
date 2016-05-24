##' Likelihood profiles for fitted strategic models
##' 
##' Calculate profile likelihood to assess convergence of a model.
##'
##' Likelihood profiling can help determine if a model fit failed to reach a
##' global maximum, which can be an issue (especially for the
##' \code{\link{ultimatum}} model).  The process of profiling is as follows: a
##' parameter selected to be profiled is fixed at certain values spaced around
##' its originally estimated value, while the log-likelihood is maximized with
##' respect to the other parameters in the model.  For models with large numbers
##' of observations or parameters, profiling may take a long time, as \eqn{p
##' \times (2s + 1)}{p * (2s + 1)} models will be fit (p: number of parameters;
##' s: number of steps).
##'
##' The function will issue a warning if a model fit in profiling has a
##' log-likelihood exceeding that of the original model.  This means the
##' original fit failed to reach a global maximum, and any inferences based on
##' the fitted model are invalid.  If this occurs, refit the model, passing the
##' \code{profile.game} output to the fitting function's \code{profile} argument
##' (as in the example below).  The new fit will use the coefficients from the
##' profile fit with the highest log-likelihood as starting values.
##'
##' The function is based loosely on \code{\link{profile.glm}} in the \pkg{MASS}
##' package.  However, that function focuses on the calculation of exact
##' confidence intervals for regression coefficients, whereas this one is for
##' diagnosing non-convergence.  Future versions of the \pkg{games} package may
##' incorporate the confidence interval functionality as well.
##' @param fitted a fitted model of class \code{game}.
##' @param which integer vector giving the indices of the parameters to be
##' profiled.  The default is to use all parameters.  Parameters that were held
##' fixed in the original fitting are ignored if selected.
##' @param steps number of steps to take (in each direction) from the original
##' value for each parameter to be profiled.
##' @param dist distance the last step should be from the original parameter
##' value  (in terms of standard errors if \code{use.se} is \code{TRUE};
##' absolute value otherwise).  Should be a numeric vector of length equal to 1
##' or \code{length(coef(fitted))}.
##' @param use.se logical: whether \code{dist} refers to standard errors
##' @param report logical: whether to print status bar (for complex models or
##' those with many parameters, profiling can be lengthy)
##' @param ... other arguments to be passed to the fitting function (see
##'ode{\link{maxLik}}).
##' @return A list of data frames, each containing the estimated coefficients
##' across the profiled values for a particular parameter.  The first column of
##' each data frame is the log-likelihood for the given fits.  The returned
##' object is of class \code{c("profile.game", "profile")}.
##' @seealso \code{\link{plot.profile.game}} for plotting profiled likelihoods
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @example inst/examples/profile.game.r
profile.game <- function(fitted, which = 1:p, steps = 5, dist = 3, use.se =
                         TRUE, report = TRUE, ...)
{
    ## get the regressors from the original model
    mf <- match(c("subset", "na.action"), names(fitted$call), 0L)
    mf <- fitted$call[c(1L, mf)]
    mf$formula <- fitted$formulas
    mf$data <- fitted$model
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    regr <- list()
    for (i in seq_len(length(fitted$formulas)[2]))
        regr[[i]] <- model.matrix(fitted$formulas, data = mf, rhs = i)

    ## other info needed to refit the original model: coefficients (for starting
    ## values), fixed parameters, outcome, link, and type.  also will retrieve
    ## the standard errors to determine values to profile at
    cf <- fitted$coefficients
    p <- length(cf)
    if (length(dim(fitted$y))) {  ## for ultimatum objects (if offer and accept
                                  ## both included)
        y <- fitted$y[, 1]
    } else {
        y <- as.numeric(fitted$y)  ## use as.numeric() because y is stored as a
                                   ## factor in a game object
    }

    fixed <- fitted$fixed
    link <- fitted$link
    type <- fitted$type
    method <- fitted$convergence$method
    if (length(dist) == 1L)
        dist <- rep(dist, length(cf))
    if (use.se) {
        se <- sqrt(diag(fitted$vcov))
    } else {
        se <- rep(1, length(cf))
    }

    ## special stuff for ultimatum models
    if (inherits(fitted, "ultimatum")) {
    	minOffer <- fitted$minOffer
        maxOffer <- fitted$maxOffer
        acc <- if (length(dim(fitted$y))) fitted$y[, 2] else NULL
        offerOnly <- fitted$outcome == "offer"
        offertol <- fitted$offertol
        names(regr) <- c("X", "Z")
    } else {
        minOffer <- maxOffer <- acc <- offerOnly <- offertol <- NULL
    }

    ## ## Special stuff for egame12cor models
    ## if (inherits(fitted, "egame12cor")) {
    ##     gridsize <- fitted$gridsize
    ## } else {
    ##     gridsize <- NULL
    ## }

    ## retrieve the appropriate log-likelihood and gradient, based on the type
    ## of model.  this relies on the class attribute of the object being ordered
    ## correctly, which will generally be the case if users don't mess with the
    ## objects created, but in the future it might be better to use something
    ## like the "family" attribute of glm objects
    logLik <- switch(class(fitted)[2],
                     egame12 = logLik12,
                     ## egame12cor = logLik12cor,
                     egame122 = logLik122,
                     egame123 = logLik123,
                     ultimatum = logLikUlt)
    logLikGrad <- switch(class(fitted)[2],
                         egame12 = logLikGrad12,
                         ## egame12cor = NULL,
                         egame122 = logLikGrad122,
                         egame123 = logLikGrad123,
                         ultimatum = logLikGradUlt)
    if (!fitted$convergence$gradient)
        logLikGrad <- NULL

    ## looping over each parameter value
    ##
    ## the code in the loop is loosely based on that of MASS:::profile.glm,
    ## written by D.M. Bates and W.N. Venables, licensed under the GPL.  the
    ## difference is that theirs uses a nice adaptive algorithm to go
    ## approximately the right distance in the specified number of steps,
    ## whereas mine is dumb and depends entirely on the supplied values of
    ## 'steps' and 'dist'
    didNotConverge <- FALSE
    ans <- list()
    k <- 0
    if (report) {
        cat("\nEstimating likelihood profiles...\n")
        pb <- txtProgressBar(min = 1, max = length(which) * (2*steps + 1))
    }
    for (i in which) {
        if (fixed[i])  # skip fixed parameters
            next

        ## calculate the new parameter values to maximize at
        thisAns <- matrix(nrow = 2*steps, ncol = length(cf) + 1)
        thisAns <- data.frame(thisAns)
        names(thisAns) <- c("logLik", names(cf))
        cfvals <- seq(cf[i] - dist[i]*se[i], cf[i] + dist[i]*se[i], length.out =
                      2*steps + 1)
        cfvals <- cfvals[-(steps + 1)]  # prevent refitting from original
                                        # coefficients (causes the
                                        # non-convergence warning to be issued
                                        # inappropriately)
        fvec <- fixed
        fvec[i] <- TRUE

        ## inner loop: refitting while fixing parameter i at the j'th element of
        ## 'cfvals'
        for (j in seq_along(cfvals)) {
            sval <- cf
            sval[i] <- cfvals[j]
            results <-
                tryCatch(maxLik(logLik = logLik, grad = logLikGrad, start =
                                sval, fixed = fvec, method = method, y = y, regr
                                = regr, link = link, type = type, acc = acc,
                                minOffer=minOffer, maxOffer = maxOffer, offerOnly = offerOnly,
                                offertol = offertol, ...),
                         error = identity)
            if (inherits(results, "error")) {
                thisAns[j, ] <- NA
            } else {
                thisAns[j, ] <- c(results$max, results$estimate)
            }

            k <- k + 1
            setTxtProgressBar(pb, k)
        }

        ## check if any of the profiled fits yield a log-likelihood higher than
        ## that of the the original fitted model
        if (isTRUE(any(thisAns[, 1] > sum(fitted$log.likelihood))))
            didNotConverge <- TRUE
        ans[[names(cf)[i]]] <- thisAns

    }
    if (report)
        cat("\n")

    if (didNotConverge)
        warning("some profiled fits have higher log-likelihood than original fit; refit the model using \"profile\" option")

    attr(ans, "original.fit") <- fitted
    class(ans) <- c("profile.game", "profile")
    return(ans)
}

##' Plot profiles of strategic model log-likelihoods
##' 
##' Plot output of \code{\link{profile.game}}.
##'
##' This method provides plots for a quick assessment of whether \code{game}
##' models have failed to converge to a global maximum.  For each parameter
##' profiled (see \code{\link{profile.game}} for details of the profiling
##' process), a spline interpolation of the log-likelihood profile is provided,
##' with an "x" marking the value at the original parameter estimate.
##'
##' Sometimes the plot will seem to indicate that the original fit did not reach
##' the global maximum, even though \code{profile.game} did not issue the
##' non-covergence warning.  This is an artifact of the interpolation, which can
##' be confirmed by re-running \code{plot.profile.game} with the argument
##' \code{show.pts = TRUE}.
##' @param x an object of class \code{profile.game}, typically created by
##' running \code{\link{profile.game}} on a fitted \code{game} model
##' @param show.pts logical: plot a point for the log-likelihood of each
##' profiled model?
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly
##' @method plot profile.game
##' @seealso \code{\link{profile.game}}
##' @export
##' @author Brenton Kenkel
plot.profile.game <- function(x, show.pts = FALSE, ...)
{
    ## this function is based largely on MASS:::plot.profile, written by
    ## D.M. Bates and W.N. Venables, licensed under GPL

    nm <- names(x)
    origcf <- coef(attr(x, "original.fit"))
    origll <- sum(attr(x, "original.fit")$log.likelihood)

    ## saves the user's original graphical parameters (which are restored after
    ## the new plot is made) and sets new ones to display all plots at once
    nr <- ceiling(sqrt(length(nm)))
    oldpar <- par(mfrow = c(nr, nr))
    on.exit(par(oldpar))

    for (nam in nm) {
        xvals <- x[[nam]][, nam]
        yvals <- x[[nam]][, "logLik"]

        ## make sure interpolation passes through original fit
        xvals <- c(xvals, origcf[nam])
        yvals <- c(yvals, origll)

        splineVals <- spline(xvals, yvals)
        ylim <- c(min(yvals, splineVals$y), max(yvals, splineVals$y))

        plot(xvals, yvals, xlab = nam, ylab = "log-likelihood", type = "n",
             ylim = ylim)

        ## plot an "x" at the main model estimate
        points(xvals[length(xvals)], yvals[length(yvals)], pch = 4)

        ## plot other points (if desired)
        if (show.pts) {
            points(xvals[-length(xvals)], yvals[-length(xvals)])
        }

        ## spline fit for the likelihood curve
        lines(splineVals$x, splineVals$y)
    }

    invisible(x)
}
