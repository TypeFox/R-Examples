##' Mode of a vector
##' 
##' Finds the modal value of a vector of any class.
##'
##' Based on the Stack Overflow answer
##' \url{http://stackoverflow.com/a/8189441/143383}
##' @param x a vector (lists and arrays will be flattened).
##' @param na.rm logical: strip \code{NA} values?
##' @return The value of \code{x} that occurs most often.  If there is a tie,
##' the one that appears first (among those tied) is chosen.
##' @export
##' @author Ken Williams (on Stack Overflow)
##' @example inst/examples/Mode.r
Mode <- function(x, na.rm = FALSE)
{
    if (na.rm)
        x <- x[!is.na(x)]
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

##
## INPUT:
## x: a data frame
## ...: expressions (see below)
##
## RETURN:
## a one-row data frame
##
## The function takes a data frame 'x' and returns a one-row data frame with the
## same variables, containing a "typical observation" profile from x.  The
## default is to take the mean of numeric variables, the median of ordered
## and binary variables, and the mode of categorical variables.
##
## These defaults can be overridden by passing expressions to "...".  For
## example, to set a variable 'z' to its 25th percentile, use:
##     makeProfile(x, z = quantile(z, 0.25))
## To set 'z' to equal 0.3 and 'w' to equal its maximum, use:
##     makeProfile(x, z = 0.3, w = max(w))
##
## This function is used in 'predProbs' and is not meant to be called by users.
## It is analogous to the 'setx' function in the Zelig package.
##
makeProfile <- function(x, ...)
{
    cl <- match.call(expand.dots = FALSE)

    ## get the first row of the data frame (and continue to store as data frame
    ## to allow for different data types, preserve factor levels, etc)
    ans <- x[1, ]

    ## loop over each variable in the data frame
    for (i in seq_len(ncol(x))) {
        xvar <- x[, i]

        ## use mean for numeric variables, median for ordered or dummy
        ## variables, mode for categorical variables
        isDummy <- all(unique(xvar) %in% c(0, 1))
        if (is.numeric(xvar) && !isDummy) {
            ans[, i] <- mean(xvar)
        } else if (is.ordered(xvar)) {
            ans[, i] <- median(xvar)
        } else if (isDummy) {
            ans[, i] <- median(xvar)
        } else {
            ans[, i] <- Mode(xvar)
        }
    }

    if ("..." %in% names(cl)) {        
        ## takes the expressions fed to "...", evaluates them within the
        ## supplied data frame, and returns them as a vector.
        ##
        ## e.g., if the call is makeProfile(x, foo = median(foo), bar =
        ## quantile(bar, .25)), this will return a vector with named elements
        ## "foo" and "bar".
        ##
        ## the eval(substitute()) business is to ensure that these aren't
        ## evaluated in the global environment rather than the data frame "x",
        ## which would most likely throw an error when "foo" and "bar" weren't
        ## found in the workspace, or possibly obtain the wrong values for
        ## these.
        dots <- eval(substitute(list(...)), x)
        dots <- unlist(lapply(dots, unname))

        toReplace <- match(names(dots), names(x), 0L)
        ans[toReplace] <- dots[toReplace != 0L]
    }

    ## ensures that choices in "..." expressed as characters (e.g., foo = "a"
    ## for a factor foo with levels a, b, c) don't wind up turning factor
    ## variables into characters
    fvars <- sapply(x, inherits, what = "factor")
    for (i in seq_len(ncol(x))) {
        if (fvars[i]) {
            ans[, i] <- factor(ans[, i], levels = levels(x[, i]))
        }
    }
    
    return(ans)
}

##
## INPUT:
## x: fitted model of class "game" containing element 'boot.matrix'
## newdata: data frame to get predicted probabilities for
## ci: width of confidence bands
## type: outcome or action
## report: whether to print status bar
##
## RETURN:
## list of bootstrapped lower and upper confidence bands ('lows' and 'highs'
## resp.) for the predicted probabilities
## 
CIfromBoot <- function(x, newdata, ci = .95, type, report = TRUE)
{
    ## the dimensions of the predicted-probability matrix varies with the type
    ## of model (3 for egame12, 4 for egame122, 2 for ultimatum, etc), so this
    ## is just to figure out the correct number
    n <- nrow(x$boot.matrix)
    forDims <- predict(x, newdata = newdata, type = type)
    ans <- vector("list", ncol(forDims))
    for (i in seq_along(ans))
        ans[[i]] <- matrix(nrow = n, ncol = nrow(forDims))

    ## calculate the predicted values for each observation using each
    ## bootstrapped coefficient vector
    if (report) {
        cat("\nCalculating confidence intervals...\n")
        pb <- txtProgressBar(min = 1, max = n)
    }
    for (i in seq_len(n)) {
        x$coefficients <- x$boot.matrix[i, ]
        xpred <- predict(x, newdata = newdata, type = type)
        for (j in seq_along(ans))
            ans[[j]][i, ] <- xpred[, j]
        if (report)
            setTxtProgressBar(pb, i)
    }
    if (report)
        cat("\n")

    ## calculate quantiles of the predicted probabilities
    q <- .5 - (ci / 2)
    lows <- lapply(ans, function(x) apply(x, 2, quantile, probs = q))
    lows <- do.call(cbind, lows)
    highs <- lapply(ans, function(x) apply(x, 2, quantile, probs = 1 - q))
    highs <- do.call(cbind, highs)

    colnames(lows) <- paste(colnames(forDims), "low", sep = ":")
    colnames(highs) <- paste(colnames(forDims), "high", sep = ":")

    return(list(lows = lows, highs = highs))
}

##' User-friendly predicted probability analysis
##' 
##' Easy generation and plotting of predicted probabilities from a fitted
##' strategic model.
##'
##' \code{predProbs} provides an easy way to analyze the estimated marginal
##' effect of an independent variable on the probability of particular outcomes,
##' using the estimates returned by a strategic model.  The procedure is
##' designed so that, for a preliminary analysis, the user can simply specify
##' the fitted model and the independent variable of interest, and quickly
##' obtain plots of predicted probabilities.  However, it is flexible enough to
##' allow for finely tuned analysis as well.
##' 
##' The procedure works by varying \code{x}, the variable of interest, across
##' its observed range (or one specified by the user in \code{xlim}) while
##' holding all other independent variables in the model fixed.  The profile
##' created by default is as follows (the same defaults as in the \code{sim}
##' function in the \pkg{Zelig} package):
##' \itemize{
##' \item numeric, non-binary variables are fixed at their means
##' \item \code{\link{ordered}} variables are fixed at their medians
##' \item all others are fixed at their modes (see \code{\link{Mode}})}
##' However, it is possible to override these defaults for any or all
##' variables.  For example, to set a variable named \code{polity} to its lower
##' quartile, call \code{predProbs} with the argument \code{polity =
##' quantile(polity, 0.25)}.  To set a factor variable to a particular level,
##' provide the name of the level as a character string (in quotes). (Also see
##' the examples below.)
##'
##' Confidence intervals for each predicted point are generated by bootstrap.
##' If \code{model} has a non-null \code{boot.matrix} element (i.e., a bootstrap
##' was performed with the model fitting), then these results are used to
##' make the confidence intervals.  Otherwise, a parametric bootstrap sample is
##' generated by sampling from a multivariate normal distribution around the
##' parameter estimates.  In this case, a warning is issued.
##'
##' For information on plotting the predicted probabilities, see
##' \code{\link{plot.predProbs}}.  The plots are made with base graphics.  If you
##' prefer to use an alternative graphics package, all the information necessary
##' to make the plots is included in the data frame returned.
##' @param model a fitted model of class \code{game}.
##' @param x character string giving the name of the variable to place "on the
##' x-axis" while all others are held constant.  Partial matches are accepted.
##' @param xlim numeric, length 2: the range that \code{x} should be varied over
##' (if \code{x} is continuous).  Defaults to the observed range of \code{x}.
##' @param n integer: the number of observations to generate (if \code{x} is
##' continuous).
##' @param ci numeric: width of the confidence interval to estimate around each
##' predicted probability.  Set to \code{0} to estimate no confidence intervals.
##' @param type whether to generate predicted values for outcomes (the default)
##' or actions
##' @param makePlots logical: whether to automatically make the default plot
##' for the returned object.  See \code{\link{plot.predProbs}}.
##' @param report logical: whether to print a status bar while obtaining the
##' confidence intervals for the predicted probabilities.
##' @param ... used to set values for variables other than \code{x} in the
##' profile of observations.  See "Details" and "Examples".
##' @return An object of class \code{predProbs}.  This is a data frame containing
##' each hypothetical observation's predicted probability, the upper and lower
##' bounds of the confidence interval, and the value of each regressor.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com}).  Code for
##' escaping special regex characters was taken from the \code{Hmisc} package's
##' function \code{escapeRegex}, written by Charles Dupont.
##' @seealso \code{\link{predict.game}} for somewhat more flexible (but fussier)
##' generation of predicted probabilities.
##' @example inst/examples/predProbs.r
predProbs <- function(model, x, xlim = c(min(x), max(x)), n = 100, ci = .95,
                      type = c("outcome", "action"),
                      makePlots = FALSE, report = TRUE, ...)
{
    type <- match.arg(type)

    ## find the variable corresponding to the named 'x' and stop if it matches
    ## none or more than one
    xc <- charmatch(x, names(model$model))
    if (is.na(xc)) {
        stop(dQuote(x),
             " does not match any variable names in the supplied model: ",
             paste(names(model$model), collapse = ", "))
    } else if (xc == 0) {
        ## The regular expression in the next line is taken from the escapeRegex
        ## function in the Hmisc package (v3.8-2, 2010-06-22), written by
        ## Charles Dupont, licensed under GPL
        xc <- paste("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", x), sep =
                    "")
        xc <- grep(xc, names(model$model), value = TRUE)
        stop(dQuote(x),
             " matches multiple variable names in the supplied model: ",
             paste(xc, collapse = ", "))
    } else {
        x <- model$model[, xc]
    }

    ## construct the set of values of 'x' to evaluate at
    if (all(unique(x) %in% c(0, 1))) {
        ## If 'x' is binary, just use 0 and 1
        xs <- c(0, 1)
    } else if (is.numeric(x)) {
        ## If 'x' is numeric, make a grid of size 'n' within 'xlim'
        xs <- seq(xlim[1], xlim[2], length.out = n)
    } else if (is.factor(x)) {
        ## If 'x' is a factor, use each of its levels
        xs <- rep(x[1], nlevels(x))
        xs[] <- levels(x)
    } else if (is.logical(x)) {
        ## If 'x' is logical, use the logical values
        xs <- c(FALSE, TRUE)
    }

    ## construct a profile of values for the variables other than 'x'; see
    ## 'makeProfile' above for how this is done
    prof <- makeProfile(model$model, ...)
    profData <- prof[rep(1, length(xs)), , drop = FALSE]
    profData[, xc] <- xs
    rownames(profData) <- seq_along(xs)

    ## get predicted probabilities from the constructed profile
    ans <- predict(model, newdata = profData, type = type)
    if (is.list(ans))
        ans <- do.call(cbind, ans)

    ## bootstrap confidence intervals for the predicted probabilities, using a
    ## Clarify-style parametric resampling if the model was not bootstrapped initially
    if (ci > 0 && is.null(model$boot.matrix)) {
        warning("Bootstrap values unavailable; using normal resampling to generate confidence intervals")
        vcv <- vcov(model)[!model$fixed, !model$fixed, drop = FALSE]
        bm <- mvrnorm(1000, mu = coef(model)[!model$fixed], Sigma = vcv)
        bbm <- matrix(NA, nrow = 1000, ncol = length(coef(model)))
        bbm[, !model$fixed] <- bm
        for (i in seq_along(model$fixed)) {
            if (model$fixed[i])
                bbm[, i] <- coef(model)[i]
        }
        model$boot.matrix <- bbm
    }
    if (ci > 0)
        CIvals <- CIfromBoot(model, newdata = profData, ci = ci, type = type)

    ## save the indices of the columns that the probabilities, confidence
    ## bounds, and variable of interest are in, and make them attributes of the
    ## output in order for 'plot.predProbs' to use them
    probcols <- 1:ncol(ans)
    if (ci > 0) {
        ans <- cbind(ans, do.call(cbind, unname(CIvals)))
        lowcols <- max(probcols) + probcols
        highcols <- max(lowcols) + probcols
    } else {
        lowcols <- highcols <- numeric(0)
    }
    xcol <- max(probcols, lowcols, highcols) + xc
    ans <- structure(cbind(ans, profData), probcols = probcols, lowcols =
                     lowcols, highcols = highcols, xcol = xcol, class =
                     c("predProbs", "data.frame"))

    if (makePlots)
        plot(ans)

    invisible(ans)
}

## -----------------------------------------------------------------------------
## NOTE:
## The next two functions are designed to mimic the behavior of the "plot"
## method for "gam" objects.  See the source of "plot.gam" and
## "plot.preplot.gam" for more.  The "gam" package is licensed under the GPL,
## and some of the code below closely matches it.
## -----------------------------------------------------------------------------

##' Plot predicted probabilities
##' 
##' Plots predicted probabilities and associated confidence bands, using the
##' data returned from a call to \code{\link{predProbs}}.
##'
##' Most \code{predProbs} objects will be associated with multiple plots: one for
##' each outcome in the estimated model.  These are the three or four terminal
##' nodes for a \code{\link{egame12}} or \code{\link{egame122}} model
##' respectively; for an \code{\link{ultimatum}} model, these are the expected
##' offer and the probability of acceptance.  By default, \code{plot.predProbs}
##' produces plots for all of them, so only the last will be visible unless the
##' graphics device is set to have multiple figures (e.g., by setting
##' \code{par(mfrow = ...)}).  The argument \code{ask} displays a menu to select
##' among the possible plots for a given object, and \code{which} allows for
##' this to be done non-interactively.
##' @param x an object of class \code{predProbs} (i.e., a data frame returned by
##' \code{\link{predProbs}}).
##' @param which optional integer specifying which plot (as numbered in the menu
##' displayed when \code{ask == TRUE}) to make.  If none is given, all available
##' plots are printed in succession.
##' @param ask logical: display interactive menu with options for which plot to
##' make?
##' @param ... further arguments to pass to the plotting function.  See
##' \code{\link{plot.default}} (when the variable on the x-axis is continuous)
##' or \code{\link{bxp}} (when it is discrete).
##' @return an object of class \code{preplot.predProbs}, invisibly.  This contains
##' the raw information used by lower-level plotting functions.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @example inst/examples/plot.predProbs.r
plot.predProbs <- function(x, which = NULL, ask = FALSE, ...)
{
    ## retrieve the columns containing the variable of interest, predicted
    ## probabilities, and confidence bands for each outcome
    probs <- x[, attr(x, "probcols"), drop = FALSE]
    if (length(attr(x, "lowcols"))) {
        lows <- x[, attr(x, "lowcols"), drop = FALSE]
        highs <- x[, attr(x, "highcols"), drop = FALSE]
    } else {
        lows <- highs <- NULL
    }
    xvar <- x[, attr(x, "xcol")]
    probnames <- names(x)[attr(x, "probcols")]

    ## make a "preplot" object, which is a list of lists containing the plot
    ## components for each outcome
    preplotObj <- vector("list", ncol(probs))
    class(preplotObj) <- "preplot.predProbs"
    for (i in 1:ncol(probs)) {
        preplotObj[[i]] <- list()
        preplotObj[[i]]$x <- xvar
        preplotObj[[i]]$y <- probs[, i]
        preplotObj[[i]]$low <- lows[, i]
        preplotObj[[i]]$high <- highs[, i]
        preplotObj[[i]]$xlab <- names(x)[attr(x, "xcol")]
        preplotObj[[i]]$ylab <- probnames[i]
        class(preplotObj[[i]]) <- "preplot.predProbs"
    }

    if (is.null(which) && ask) {
        ## run an interactive menu for choosing which plot to view
        tmenu <- c(paste("plot:", probnames), "plot all terms")
        pick <- 1
        while (pick > 0 && pick <= length(tmenu)) {
            pick <- menu(tmenu, title =
                         "Make a plot selection (or 0 to exit):\n")
            if (pick > 0 && pick < length(tmenu)) {
                plot.preplot.predProbs(preplotObj[[pick]], ...)
            } else if (pick == length(tmenu)) {
                plot.preplot.predProbs(preplotObj, ...)
            }
        }
    } else if (!is.null(which)) {
        ## plot the requested outcomes
        plot.preplot.predProbs(preplotObj[[which]], ...)
    } else {
        ## plot all outcomes
        plot.preplot.predProbs(preplotObj, ...)
    }

    invisible(preplotObj)
}

##' @method plot preplot.predProbs
##' @export
plot.preplot.predProbs <- function(x, xlab = x$xlab, ylab = x$ylab,
                                   ylim = c(min(x$y, x$low), max(x$y, x$high)),
                                   type = "l", lty.ci = 2, ...)
{
    listof <- inherits(x[[1]], "preplot.predProbs")
    cl <- match.call()
    
    if (listof) {
        ## if 'x' contains a list of preplot objects, run the function on each
        ## of its components (producing length(x) plots in sequence)
        for (i in seq_along(x)) {
            icall <- cl
            icall$x <- x[[i]]
            eval(icall, parent.frame())
        }
    } else if (is.factor(x$x)) {
        ## if the variable of interest is a factor, then we need to make
        ## boxplots "manually" via the bxp function (see ?boxplot and ?bxp)
        boxStats <- list()
        boxStats$stats <- matrix(x$y, nrow = 5, ncol = length(x$y), byrow =
                                 TRUE)
        if (!is.null(x$low)) {
            boxStats$stats[1, ] <- x$low
            boxStats$stats[5, ] <- x$high
        }
        boxStats$n <- rep(1, length(x$x))
        boxStats$conf <- boxStats$stats[c(1, 5), ]
        boxStats$out <- numeric(0)
        boxStats$group <- numeric(0)
        boxStats$names <- as.character(x$x)

        bxp(z = boxStats, xlab = xlab, ylab = ylab, ...)
    } else {
        ## if the variable of interest is numeric, just make a typical plot
        plot(x$x, x$y, type = type, xlab = xlab, ylab = ylab, ylim = ylim, ...)
        if (!is.null(x$low)) {
            lines(x$x, x$low, lty = lty.ci)
            lines(x$x, x$high, lty = lty.ci)
        }
    }
}
