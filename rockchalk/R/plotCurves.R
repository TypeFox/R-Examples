##' Assists creation of predicted value curves for regression models.
##'
##'
##' Creates a predicted value plot that includes a separate predicted
##' value line for each value of a focal variable. The x axis variable
##' is specified by the \code{plotx} argument. As of rockchalk 1.7.x,
##' the moderator argument, modx, is optional. Think of this a new
##' version of R's \code{termplot}, but it allows for
##' interactions. And it handles some nonlinear transformations more
##' gracefully than termplot.
##'
##' This is similar to \code{plotSlopes}, but it accepts regressions
##' in which there are transformed variables, such as "log(x1)".
##' It creates a plot of the predicted dependent
##' variable against one of the numeric predictors, \code{plotx}. It
##' draws a predicted value line for each value of \code{modx}, a
##' moderator variable. The moderator may be a numeric or categorical
##' moderator variable.
##'
##' The user may designate which particular values of the moderator
##' are used for calculating the predicted value lines.  That is,
##' \code{modxVals = c( 12,22,37)} would draw lines for values 12, 22,
##' and 37 of the moderator. User may instead supply a character
##' string to choose one of the built in algorithms. The default
##' algorithm is "quantile", which will select \code{n} values that
##' are evenly spaced along the \code{modx} axis. The algorithm
##' "std.dev" will select the mean of \code{modx} (m) and then it will
##' select values that step away from the mean in standard deviation
##' sd units. For example, if \code{n = 3}, the focal
##' values will \code{m, m - sd, am + sd}.
##'
##' @param model Required. Fitted regression object. Must have a
##' predict method
##' @param plotx Required. String with name of predictor for the x axis
##' @param plotxRange Optional. If not specified, the observed
##' range of plotx will be used to determine the axis range.
##' @param nx Number of values of plotx at which to calculate the predicted
##' value.  Default = 40.
##' @param modx Optional. String for moderator variable name. May be
##' either numeric or factor.
##' @param n Optional.  Number of focal values of \code{modx}, used by
##' algorithms specified by modxVals; will be ignored if modxVals
##' supplies a vector of focal values.
##' @param modxVals Optional. A vector of focal values for which
##' predicted values are to be plotted. May also be a character string
##' to select an algorithm ("quantile","std.dev." or "table"), or a
##' user-supplied function to select focal values (a new method
##' similar to \code{getFocal}). If modx is a factor, currently, the
##' only available algorithm is "table" (see \code{getFocal.factor}.
##' @param interval Optional. Intervals provided by the
##' \code{predict.lm} may be supplied, either "conf" (95% confidence
##' interval for the estimated conditional mean) or "pred" (95%
##' interval for observed values of y given the rest of the model).
##' @param plotPoints Optional. TRUE or FALSE: Should the plot include
##' the scatterplot points along with the lines.
##' @param plotLegend Optional. TRUE or FALSE: Should the default
##' legend be included?
##' @param legendTitle Optional. You'll get an automatically generated
##' title, such as "Moderator: modx", but if you don't like that,
##' specify your own string here.
##' @param legendPct Default = TRUE. Variable labels print with sample percentages.
##' @param col Optional.  A color vector to differentiate the moderator
##' values in the plot. If not specified, the R's builtin palette()
##' will be used. User may supply a vector of valid color names,
##' either explicitly c("pink","black", "gray70") or implicitly,
##' rainbow(10) or gray.colors(5). Color names will be recycled if there
##' are more focal values of \code{modx} than colors provided.
##' @param envir environment to search for variables.
##' @param llwd Optional. Line widths for predicted values. Can be
##' single value or a vector, which will be recycled as necessary.
##' @param opacity Optional, default = 100. A number between 1 and
##' 255. 1 means "transparent" or invisible, 255 means very dark.  the
##' darkness of confidence interval regions
##' @param ... further arguments that are passed to plot or
##' predict. The arguments that are monitored to be sent to predict
##' are c("type", "se.fit", "dispersion", "interval", "level",
##' "terms", "na.action").
##' @export
##' @import car
##' @return A plot is created as a side effect, a list is returned
##' including 1) the call, 2) a newdata object that includes
##' information on the curves that were plotted, 3) a vector modxVals,
##' the values for which curves were drawn.
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @example  inst/examples/plotCurves-ex.R
plotCurves <-
    function (model, plotx, nx = 40,  modx,  plotxRange = NULL, n, modxVals = NULL,
              interval = c("none", "confidence", "prediction"),
              plotPoints = TRUE, plotLegend = TRUE, legendTitle = NULL,
              legendPct = TRUE, 
              col = NULL, llwd = 2, opacity = 100,
              envir = environment(formula(model)), ...)
{
    if (missing(model))
        stop("plotCurves requires a fitted regression model.")
    if (missing(plotx))
        stop("plotCurves requires the name of the variable to be drawn on the x axis")

    cl <- match.call()
    interval <- match.arg(interval)
    mf <- model.frame(model)
    emf <- model.data(model)

    zzz <- as.character(substitute(plotx))
    plotx <- zzz[1L]
    
    if (!missing(modx)) {
        zzz <- as.character(substitute(modx))
        modx <- zzz[1L]
    }
        
    plotxVar <- emf[ , plotx]
    if (!is.numeric(plotxVar))
        stop(paste("plotCurves: The variable", plotx, "should be a numeric variable"))

    ## Gamble on row names to select which cases are nonmissing
    depVar <- model.response(mf)[row.names(emf)]

    ylab <- names(mf)[1]  ## returns transformed DV
    plotxRange <- if(is.null(plotxRange)) range(plotxVar, na.rm = TRUE) else plotxRange
    plotxVals <- plotSeq(plotxRange, length.out = nx)

    ##Bug noticed 2013-09-22
    ## focalVals for a nonlinear model needs to be long set,
    ## not just end points as I imagined before. See 3 lines below
    if (missing(modx) || is.null(modx)) {
        modxVar <- rep(1, nobs(model))
        if (interval == "none") {
            ## focalVals <- list(plotxRange)
            focalVals <- list(plotxVals)
        } else {
            focalVals <- list(plotxVals)
        }
        names(focalVals) <- c(plotx)
        modx <- NULL
        modxVals <- 1
    } else {
        modxVar <- emf[ , modx]
        if (is.factor(modxVar)) { ## modxVar is a factor
            n <- ifelse(missing(n), nlevels(modxVar), n)
            modxVals <- getFocal(modxVar, xvals = modxVals, n, pct = legendPct)
        } else {
            n <- ifelse(missing(n), 3, n)
            modxVals <- getFocal(modxVar, xvals = modxVals, n, pct = legendPct)
        }

        focalVals <- list(modxVals, plotxVals)
        names(focalVals) <- c(modx, plotx)
    }

    newdf <- newdata(model, predVals = focalVals, emf = emf)

    dotargs <- list(...)
    dotnames <- names(dotargs)
    ## scan dotargs for predict keywords. Remove from dotargs
    ## the ones we only want going to predict. Leave
    ## others.

    level <- if (!is.null(dotargs[["level"]])) dotargs[["level"]] else 0.95

    parms <- list(model, newdata = newdf, type = "response" , interval = interval)
    ## type requires special handling because it may go to predict or plot
    if (!is.null(dotargs[["type"]])) {
        if (dotargs[["type"]] %in% c("response", "link", "none")) {
            parms[["type"]] <- dotargs[["type"]]
            dotargs[["type"]] <- NULL
        }
    }
    predArgs <- list()
    validForPredict <- c("se.fit", "dispersion", "level", "terms", "na.action")
    dotsForPredict <- dotnames[dotnames %in% validForPredict]

    if (length(dotsForPredict) > 0){
        parms <- modifyList(parms, dotargs[dotsForPredict])
        dotargs[[dotsForPredict]] <- NULL
    }
    np <- do.call("predictCI", parms)
    newdf <- cbind(newdf, np$fit)

    if ((!is.null(parms[["se.fit"]])) && (parms[["se.fit"]] == TRUE)) newdf <- cbind(newdf, np$se.fit)

    if ("ylim" %in% dotargs) {
        plotyRange <- dotargs$ylim
        dotargs$ylim <- NULL
    } else if (is.logical(depVar) || (is.factor(depVar) && length(levels(depVar)) == 2)) {
        plotyRange <- c(0, 1.2)
    } else if (is.numeric(depVar)) {
        plotyRange <- magRange(depVar, mult = c(1, 1.2))
    } else {
        stop("plotCurves: The dependent variable is neither numeric nor logical. I don't know what you want me to do. Please be patient, I'll figure it out")
    }
    
    parms <- list(newdf = newdf, olddf = data.frame(modxVar, plotxVar, depVar),
                  plotx = plotx, modx = modx, modxVals = modxVals,
                  interval = interval, level = level, plotPoints = plotPoints,
                  plotLegend = plotLegend, legendTitle = legendTitle,
                  col = col,  opacity = opacity, xlim = plotxRange,
                  ylab = ylab, ylim = plotyRange, llwd = llwd)
    
    parms <- modifyList(parms, dotargs)
    plotArgs <- do.call("plotFancy", parms)

    z <- list(call = cl, newdata = newdf, modxVals = modxVals,
              col = plotArgs$col, lty = plotArgs$lty)

    class(z) <- "rockchalk"

    invisible(z)
}
