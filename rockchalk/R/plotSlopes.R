##' Generic function for plotting regressions and interaction effects
##'
##' This is a generic function for plotting regression
##' objects. So far, there is an implementation for \code{lm()} objects.
##' This allows interaction effects, but not nonlinearities like log(x1).
##' For that, please see \code{plotCurves}.
##'
##' @param model Required. A fitted Regression
##' @param plotx Required. Name of one predictor from the fitted
##' model to be plotted on horizontal axis
##' @param ... Additional arguments passed to methods. Often includes
##' arguments that are passed to plot. Any arguments that customize
##' plot output, such as lwd, cex, and so forth, may be
##' supplied. These arguments intended for the predict method will be
##' used: c("type", "se.fit", "interval", "level", "dispersion",
##' "terms", "na.action")
##' @export plotSlopes
##' @importFrom graphics abline arrows legend lines mtext plot points polygon text
##' @importFrom methods as
##' @importFrom stats AIC addmargins as.formula coef cor
##'             df.residual drop1 family fitted formula lm
##'             model.frame model.matrix model.response na.omit na.pass
##'             napredict nobs pchisq pf predict
##'             predict.lm pt qnorm qt quantile resid
##'             rnorm sd setNames terms var vcov
##' @importFrom utils browseURL methods modifyList 
##' @rdname plotSlopes
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @seealso \code{\link[rockchalk]{testSlopes}} \code{\link[rockchalk]{plotCurves}}
##' @return Creates a plot and an output object that summarizes it.
plotSlopes <- function(model, plotx, ...) UseMethod("plotSlopes")

##' Plot predicted values for focal values of a moderator variable.
##'
##' This is a "simple slope" plotter for linear regression objects
##' that are created by \code{lm()}.  The function \code{plotCurves()}
##' can handle nonlinear predictive equations and generalized linear
##' models. The term "simple slopes" was coined by psychologists
##' (Aiken and West, 1991; Cohen, et al 2002) for analysis of
##' interaction effects for particular values of a moderating
##' variable. The moderating variable may be continuous or
##' categorical, lines will be plotted for focal values of that
##' variable.
##'
##' This function works well with lm models in which the predictor
##' formula includes interactions, but it does not work well with
##' nonlinear predictors (log(x) and poly(x)).  For that, please use
##' \code{plotCurves}. plotSlopes is needed only when one wants to
##' create an output object that can be used as input for
##' \code{testSlopes()}.
##'
##' The argument \code{plotx} is the name of the horizontal plotting
##' variable; it must be numeric.  The argument \code{modx} is the
##' moderator variable. It may be either a numeric or a factor
##' variable. As of version 1.7, the modx argument may be omitted. A
##' single predicted value line will be drawn. That version also
##' introduced the arguments interval and n.
##'
##' There are many ways to specify focal values using the arguments
##' \code{modxVals} and \code{n}. This changed in rockchalk-1.7.0.  If
##' \code{modxVals} is omitted, a default algorithm will be used,
##' selecting \code{n} values for plotting. \code{modxVals} may be a
##' vector of values (for a numeric moderator) or levels (for a
##' factor).  If modxVals is a vector of values, then the argument
##' \code{n} is ignored.  However, if modxVals is one of the name of
##' one of the algorithms, "table", "quantile", or "std.dev.", then
##' the argument \code{n} sets number of focal values to be selected.
##' For numeric \code{modx}, n defaults to 3, but for factors
##' \code{modx} will be the number of observed values of
##' \code{modx}. If modxVals is omitted, the defaults will be used
##' ("table" for factors, "quantile" for numeric variables).
##'
##' For the predictors besides \code{modx} and \code{plotx} (the ones
##' that are not explicitly included in the plot), predicted values
##' are calculated with variables set to the mean and mode, for numeric
##' or factor variables (respectively). Those values can be reviewed
##' in the newdata object that is created as a part of the output from
##' this function
##'
##' @param plotxRange Optional. If not specified, the observed
##' range of plotx will be used to determine the axis range.
##' @param modx Optional. String for moderator variable name. May be
##' either numeric or factor. If omitted, a single predicted value line
##' will be drawn.
##' @param n Optional. Number of focal values of
##' \code{modx}, used by algorithms specified by modxVals; will be
##' ignored if modxVals supplies a vector of focal values.
##' @param modxVals Optional. Focal values of \code{modx} for which
##' lines are desired. May be a vector of values or the name of an
##' algorithm, "quantile", "std.dev.", or "table".
##' @param interval Optional. Intervals provided by the
##' \code{predict.lm} may be supplied, either "confidence" (confidence
##' interval for the estimated conditional mean) or "prediction"
##' (interval for observed values of y given the rest of the model).
##' The level can be specified as an argument (which goes into ...
##' and then to the predict method)
##' @param plotPoints Optional. TRUE or FALSE: Should the plot include
##' the scatterplot points along with the lines.
##' @param plotLegend Optional. TRUE or FALSE: Include a default
##' @param legendTitle Optional. You'll get an automatically generated
##' title, such as "Moderator: modx", but if you don't like that,
##' specify your own string here.  legend. Set to FALSE if user wants
##' to customize a legend after the plot has been drawn.
##' @param legendPct Default = TRUE. Variable labels print with sample percentages.
##' @param col Optional. A color vector for predicted value lines (and
##' intervals if requested). If not specified, the R's builtin palette()
##' will be used. User may supply a vector of valid color names,
##' either explicitly c("pink","black", "gray70") or implicitly,
##' rainbow(10) or gray.colors(5). Color names will be recycled if there
##' are more focal values of \code{modx} than colors provided.
##' @param llwd Optional, default = 2. Line widths for predicted values. Can be
##' single value or a vector, which will be recycled as necessary.
##' @param opacity Optional, default = 100. A number between 1 and 255.
##' 1 means "transparent" or invisible, 255 means very dark.
##' Determines the darkness of confidence interval regions
##' @export
##' @method plotSlopes lm
##' @rdname plotSlopes
##' @import car
##' @return The return object includes the "newdata" object that was
##' used to create the plot, along with the "modxVals" vector, the
##' values of the moderator for which lines were drawn, and the color
##' vector. It also includes the call that generated the plot.
##' @references
##' Aiken, L. S. and West, S.G. (1991). Multiple Regression:
##' Testing and Interpreting Interactions. Newbury Park, Calif: Sage Publications.
##'
##' Cohen, J., Cohen, P., West, S. G., and Aiken, L. S. (2002).
##' Applied Multiple Regression/Correlation Analysis for the Behavioral
##' Sciences (Third.). Routledge Academic.
##' @example inst/examples/plotSlopes-ex.R
plotSlopes.lm <-
    function (model, plotx, modx, n = 3, modxVals = NULL ,
              plotxRange = NULL, interval = c("none", "confidence", "prediction"),
              plotPoints = TRUE, plotLegend = TRUE, legendTitle = NULL, legendPct = TRUE, col = NULL,
              llwd = 2, opacity = 100, ...)
{
    if (missing(model))
        stop("plotSlopes requires a fitted regression model.")
    if (missing(plotx))
        stop("plotSlopes requires the name of the variable to be drawn on the x axis")

    cl <- match.call()
    mm <- model.matrix(model)
    interval <- match.arg(interval)

    depVar <- model.response(model.frame(model))
    
    zzz <- as.character(substitute(plotx))
    plotx <- zzz[1L]
   
    if (!missing(modx)) {
        zzz <- as.character(substitute(modx))
        modx <- zzz[1L]
    }
    
    plotxVar <- model$model[, plotx]

    if (!is.numeric(plotxVar))
        stop(paste("plotSlopes: The variable", plotx, "should be a numeric variable"))
    ylab <- colnames(model$model)[1]

    plotxRange <- if(is.null(plotxRange)) range(mm[, plotx], na.rm = TRUE) else plotxRange
    plotxVals <- plotSeq(plotxRange, length.out = 40)

    ## Create focalVals object, needed by newdata
    if (missing(modx) || is.null(modx)) {
        modxVar <- rep(1, nobs(model))
        if (interval == "none") {
            focalVals <- list(plotxRange)
        } else {
            focalVals <- list(plotxVals)
        }
        names(focalVals) <- c(plotx)
        modxVals <- 1
        modx <- NULL
    } else {
        modxVar <- model$model[, modx]
        if (is.factor(modxVar)) { ## modxVar is a factor
            n <- ifelse(missing(n), nlevels(modxVar), n)
            modxVals <- getFocal(modxVar, xvals = modxVals, n, pct = legendPct)
        } else {
            n <- ifelse(missing(n), 3, n)
            modxVals <- getFocal(modxVar, xvals = modxVals, n, pct = legendPct)
        }

        ## if no interval plot requested, we only need 2 points from plotx
        ## to plot lines
        if (interval == "none") {
            focalVals <- list(modxVals, plotxRange)
        } else {
            focalVals <- list(modxVals, plotxVals)
        }
        names(focalVals) <- c(modx, plotx)
    }

    newdf <- newdata(model, predVals = focalVals)

    dotargs <- list(...)
    dotnames <- names(dotargs)
    level <- if (!is.null(dotargs[["level"]])) dotargs[["level"]] else 0.95
    
    ## scan dotargs for predict keywords. Remove from dotargs
    ## the ones we only want going to predict. Leave
    ## others.
    parms <- list(model, newdata = newdf, type = "response" , interval = interval)
    predArgs <- list()

    ## type requires special handling because it may go to predict or plot
    if (!is.null(dotargs[["type"]])) {
        if (dotargs[["type"]] %in% c("response", "link", "none")) {
            parms[["type"]] <- dotargs[["type"]]
            dotargs[["type"]] <- NULL
        }
    }

    validForPredict <- c("se.fit", "dispersion", "terms", "na.action",
                         "level", "pred.var", "weights")

    dotsForPredict <- dotnames[dotnames %in% validForPredict]

    if (length(dotsForPredict) > 0) {
        parms <- modifyList(parms, dotargs[dotsForPredict])
        dotargs[[dotsForPredict]] <- NULL
    }

    np <- do.call("predictCI", parms)
    newdf <- cbind(newdf, np$fit)
    if ((!is.null(parms[["se.fit"]])) && (parms[["se.fit"]] == TRUE)) newdf <- cbind(newdf, np$se.fit)

    plotyRange <- if(is.numeric(depVar)){
        magRange(depVar, mult = c(1, 1.2))
    } else {
        stop(paste("plotSlopes: I've not decided yet what
                   should be done when the dependent variable is not numeric.",
                   "Please be patient, I'll figure it out"))
    }

    parms <- list(newdf = newdf, olddf = data.frame(modxVar, plotxVar, depVar),
                  plotx = plotx, modx = modx, modxVals = modxVals,
                  interval = interval, level = level, plotPoints = plotPoints,
                  plotLegend = plotLegend, legendTitle = legendTitle,
                  col = col, opacity = opacity, xlim = plotxRange,
                  ylim = plotyRange, ylab = ylab, llwd = llwd)

    parms <- modifyList(parms, dotargs)
    plotArgs <- do.call("plotFancy", parms)

    z <- list(call = cl, newdata = newdf, modxVals = modxVals, col = plotArgs$col, lty = plotArgs$lty)
    class(z) <- c("plotSlopes", "rockchalk")

    invisible(z)
}
NULL

##' Regression plots with predicted value lines, confidence intervals, color coded interactions
##'
##' This is the back-end for the functions plotSlopes and plotCurves. Don't use it directly.
##'
##' @param newdf The new data frame with predictors and fit, lwr, upr variables
##' @param olddf A data frame with variables(modxVar, plotxVar, depVar)
##' @param plotx Character string for name of variable on horizontal axis
##' @param modx  Character string for name of moderator variable.
##' @param modxVals Values of moderator for which lines are desired
##' @param interval TRUE or FALSE: want confidence intervals?
##' @param plotPoints TRUE or FALSE: want to see observed values in plot?
##' @param plotLegend TRUE or FALSE: draw defaut legend
##' @param legendTitle Optional. You'll get an automatically generated title, such as "Moderator: modx",
##' but if you don't like that, specify your own string here.
##' @param col requested color scheme for lines and points. One per value of modxVals.
##' @param llwd requested line width, will re-cycle.
##' @param opacity Value in 0, 255 for darkness of interval shading
##' @param ... Other arguments passed to plot function.
##' @return col, lty, and lwd information
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##'
plotFancy <-
    function(newdf, olddf, plotx, modx, modxVals, interval, plotPoints,
             plotLegend, legendTitle, col = NULL, llwd = 2, opacity, ...)
{
    dotargs <- list(...)
    
    ## Damn. Need ylab from dotargs explicitly
    ylab <- "Dependent Variable"
    if (!is.null(dotargs[["ylab"]])) ylab <- dotargs[["ylab"]]
    ## Also need level
    if (!is.null(dotargs[["level"]])) {
        level <- dotargs[["level"]]
        dotargs[["level"]] <- NULL
    } else {
        level <- 0.95
    }
    modxVar <- olddf$modxVar
    plotxVar <- olddf$plotxVar
    depVar <- olddf$depVar

    ## Now begin the plotting work.
    if (missing(modx) || is.null(modx)) {
        lmx <- 1
    } else {
        lmx <- length(modxVals)
    }

    ## if modx is a factor's name, we want to use all the levels
    ## to set the color scheme, even if some are not used in this
    ## particular plot.
    if (is.factor(modxVar)) {
        modxLevels <- levels(modxVar)
    } else {
        modxLevels <- modxVals
        if (is.null(names(modxVals))) names(modxVals) <- modxVals
    }
    ## Deal w colors
    if (is.null(col)) {
        if (is.factor(modxVar)) {
            col <- seq_along(modxLevels)
            names(col) <- modxLevels
        }  else {
            col <- 1:lmx
            names(col) <- names(modxVals)
        }
    } else {
        if (length(col) == lmx & is.null(names(col))) {
            if (is.factor(modxVar)) names(col) <- modxLevels
            else names(col) <- names(modxVals)
        } else if (length(col) < lmx) {
            stop("plotFancy: wrong number of colors. please fix the col argument.")
        } else if (length(col) < length(modxLevels)) {
            if (is.null(names(col))) {
                names(col) <- names(modxLevels[1:length(col)])
            }
            col <- rep(col, length.out = length(modxLevels))
        }
    }

    ## Deal w line widths
    if (length(llwd) < length(col)) {
        llwd <- rep(llwd, length.out = length(col))
    }
    names(llwd) <- names(col)

    lty <- seq_along(col)
    if (!is.null(dotargs[["lty"]])) {
        lty <-  rep(dotargs[["lty"]], length.out = lmx)
        dotargs[["lty"]] <- NULL ## erase
    }
    ## lty <- seq_along(col)
    ## if (is.factor(modxVar)) {
    ##     lty <- seq_along(modxLevels)
    ## } else {
    ##     lty <- seq_along(modxVals)
    ## }
    names(lty) <- names(col)

    parms <- list(x = plotxVar, y = depVar, xlab = plotx, ylab = ylab,
                  type = "n")
    parms <- modifyList(parms, dotargs)

    do.call("plot", parms)

    ## iCol: rgb color matrix. Why does rgb insist the columns be
    iCol <- col2rgb(col)
    ### bCol: border color

    bCol <- mapply(rgb, red = iCol[1,], green = iCol[2,],
                   blue = iCol[3,], alpha = opacity, maxColorValue = 255)
    ### sCol: shade color
    sCol <-  mapply(rgb, red = iCol[1,], green = iCol[2,],
                    blue = iCol[3,], alpha = opacity/3, maxColorValue = 255)


    if (interval != "none") {
        for (j in modxVals) {
            k <- match(j, modxVals)   ##integer index
            if (is.factor(modxVar)) i <- j  ## level names
            else i <- k  ## i integer

            if (missing(modx) || is.null(modx)) {
                pdat <- newdf
            } else {
                pdat <- newdf[newdf[ , modx] %in% j, ]
            }
            parms <- list(x = c(pdat[, plotx], pdat[NROW(pdat):1 , plotx]),
                          y = c(pdat$lwr, pdat$upr[NROW(pdat):1]), lty = lty[i])
            parms <- modifyList(parms, dotargs)
            parms <- modifyList(parms, list(border = bCol[i],
                                            col = sCol[i], lwd = 0.3* llwd[k]))
            do.call("polygon", parms)
        }
    }

    for (j in modxVals) {
        if (is.factor(modxVar)) i <- j  ## level names
        else i <- match(j, modxVals)   ##integer index
        if (missing(modx) || is.null(modx)) {
            pdat <- newdf
        } else {
            pdat <- newdf[newdf[ , modx] %in% j, ]
        }
        parms <- list(x = pdat[, plotx], y = pdat$fit, lty = lty[i])
        parms <- modifyList(parms, dotargs)
        parms <- modifyList(parms, list(col = col[i], lwd = llwd[i]))
        do.call("lines", parms)
    }


    if (plotPoints) {
        parms <- list(xlab = plotx, ylab = ylab,
                      cex = 0.6, lwd = 0.75)
        if (is.factor(modxVar)) {
            parms[["col"]] <- col[as.vector(modxVar[modxVar %in% modxVals])]
            parms[["x"]] <- plotxVar[modxVar %in% modxVals]
            parms[["y"]] <- depVar[modxVar %in% modxVals]
        } else {
            parms[["col"]] <- 1
            parms[["x"]] <- plotxVar
            parms[["y"]] <- depVar
        }
        parms <- modifyList(parms, dotargs)
        do.call("points", parms)
    }

    if (plotLegend) {
        if (is.factor(modxVar)){ ## level names
            col <- col[as.vector(modxVals)]
            lty <- lty[as.vector(modxVals)]
            llwd <- llwd[as.vector(modxVals)]
        } else {
            col <- col[names(modxVals)]
            lty <- lty[names(modxVals)]
            llwd <- llwd[names(modxVals)]
        }
        if (missing(modx) || is.null(modx)) {
            titl <-  if(missing(legendTitle)) "Regression analysis"
            legnd <- c("Predicted values")
            if (interval != "none") {
                legnd[2] <- paste(level, interval, "interval")
                col <- c(col, 0)
                lty <- c(lty, 0)
                llwd <- c(llwd, 0)
            }
        } else if (is.null(names(modxVals))) {
            titl <- if(missing(legendTitle)) paste("Moderator:", modx) else legendTitle
            legnd <- paste(modxVals, sep = "")
        } else {
            titl <- if(missing(legendTitle)) paste("Moderator:", modx) else legendTitle
            legnd <- paste(names(modxVals), sep = "")
        }
        legend("topleft", legend = legnd, lty = lty, col = col,
               lwd = llwd, bg = "white", title = titl)
    }
    invisible(list(col = col, lty = lty, lwd = llwd))
}
