##' Draw a 3-D regression plot for two predictors from any linear or nonlinear lm or glm object
##'
##' This allows user to fit a regression model with many variables and
##' then plot 2 of its predictors and the output plane for those
##' predictors with other variables set at mean or mode (numeric or
##' factor).  This is a front-end (wrapper) for R's persp function.
##' Persp does all of the hard work, this function reorganizes the
##' information for the user in a more readily understood way.  It
##' intended as a convenience for students (or others) who do not
##' want to fight their way through the details needed to use persp to
##' plot a regression plane.  The fitted model can have any number of
##' input variables, this will display only two of them. And, at least
##' for the moment, I insist these predictors must be numeric
##' variables. They can be transformed in any of the usual ways, such
##' as poly, log, and so forth.
##'
##' Besides a fitted model object, plotPlane requires two additional
##' arguments, plotx1 and plotx2. These are the names of the plotting
##' variables. Please note, that if the term in the regression is
##' something like poly(fish,2) or log(fish), then the argument to
##' plotx1 should be the quoted name of the variable "fish".
##' plotPlane will handle the work of re-organizing the information so
##' that R's predict functions can generate the desired information.
##' This might be thought of as a 3D version of "termplot", with a
##' significant exception. The calculation of predicted values depends
##' on predictors besides plotx1 and plotx2 in a different ways. The
##' sample averages are used for numeric variables, but for factors
##' the modal value is used.
##'
##' This function creates an empty 3D drawing and then fills in the
##' pieces. It uses the R functions \code{lines}, \code{points}, and
##' \code{arrows}. To allow customization, several parameters are
##' introduced for the users to choose colors and such. These options
##' are prefixed by "l" for the lines that draw the plane, "p" for the
##' points, and "a" for the arrows. Of course, if plotPoints=FALSE or
##' drawArrows=FALSE, then these options are irrelevant.
##'
##' @param model an lm or glm fitted model object
##' @param plotx1 name of one variable to be used on the x1 axis
##' @param plotx2 name of one variable to be used on the x2 axis
##' @param drawArrows draw red arrows from prediction plane toward observed values TRUE or FALSE
##' @param plotPoints Should the plot include scatter of observed scores?
##' @param npp number of points at which to calculate prediction
##' @param x1lab optional label
##' @param x2lab optional label
##' @param ylab optional label
##' @param x1lim optional lower and upper bounds for x1, as vector like c(0,1)
##' @param x2lim optional lower and upper bounds for x2, as vector like c(0,1)
##' @param x1floor Default=5. Number of "floor" lines to be drawn for variable x1
##' @param x2floor Default=5. Number of "floor" lines to be drawn for variable x2
##' @param pch plot character, passed on to the "points" function
##' @param pcol color for points, col passed to "points" function
##' @param plwd line width, lwd passed to "points" function
##' @param pcex character expansion, cex passed to "points" function
##' @param llwd line width, lwd passed to the "lines" function
##' @param lcol line color, col passed to the "lines" function
##' @param llty line type, lty passed to the "lines" function
##' @param acol color for arrows, col passed to "arrows" function
##' @param alty arrow line type, lty passed to the "arrows" function
##' @param alwd arrow line width, lwd passed to the "arrows" function
##' @param alength arrow head length, length passed to "arrows" function
##' @param linesFrom object with information about "highlight" lines to be added to the 3d plane (output from plotCurves or plotSlopes)
##' @param lflwd line widths for linesFrom highlight lines
##' @param envir environment from whence to grab data
##' @param ... additional parameters that will go to persp
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @rdname plotPlane
##' @export plotPlane
##' @importFrom graphics lines
##' @seealso \code{\link[graphics]{persp}}, \code{\link[scatterplot3d]{scatterplot3d}}, \code{\link[HH]{regr2.plot}}
##' @example inst/examples/plotPlane-ex.R

plotPlane <-
function(model = NULL,  plotx1 = NULL, plotx2 = NULL, drawArrows = FALSE,
         plotPoints = TRUE, npp = 20, x1lab, x2lab, ylab, x1lim, x2lim,
         x1floor = 5, x2floor = 5,  pch = 1, pcol = "blue", plwd = 0.5,
         pcex = 1, llwd = 0.3, lcol = 1, llty = 1, acol = "red", alty = 4,
         alwd = 0.3, alength = 0.1, linesFrom,
         lflwd = 3, envir = environment(formula(model)),  ...){
    UseMethod("plotPlane")
}


##' @return The main point is the plot that is drawn, but for record
##' keeping the return object is a list including 1) res: the
##' transformation matrix that was created by persp 2) the call that
##' was issued, 3) x1seq, the "plot sequence" for the x1 dimension, 4)
##' x2seq, the "plot sequence" for the x2 dimension, 5) zplane, the
##' values of the plane corresponding to locations x1seq and x2seq.
##' @rdname plotPlane
##' @method plotPlane default
##' @export

plotPlane.default <-
function(model = NULL, plotx1 = NULL, plotx2 = NULL, drawArrows = FALSE,
         plotPoints = TRUE, npp = 20, x1lab, x2lab, ylab, x1lim, x2lim, 
         x1floor = 5, x2floor = 5, pch = 1, pcol = "blue", plwd = 0.5,
         pcex = 1, llwd = 0.3, lcol = 1, llty = 1, acol = "red", alty = 4,
         alwd = 0.3, alength = 0.1, linesFrom, lflwd = 3,
         envir = environment(formula(model)), ...){


    if (is.null(model))
        stop("plotPlane requires a fitted regression model.")
    if (is.null(plotx1) | is.null(plotx2))
        stop("plotPlane requires the name of the variable to be drawn on the x axis")
    if (plotx1 == plotx2) stop("the two plotting variables should not be the same")
    cl <- match.call()
    mf <- model.frame(model) ##y first, no intercept

    ## The dependent variable data column
    y <- model.response(mf)
    ## Need a dataframe of raw data
    emf <- model.data(model)
    varnames <- attr(emf, "varNamesRHS")

    if (!plotx1 %in% varnames) stop(paste("plotx1 variable,", plotx1, ", is not in fitted model"))
    if (!plotx2 %in% varnames) stop(paste("plotx2 variable,", plotx2, ", is not in fitted model"))
                                    
    if (!missing(x1lim))
        emf <- subset(emf, emf[ , plotx1] >= x1lim[1] & emf[ , plotx1] <= x1lim[2])
    if (!missing(x2lim))
        emf <- subset(emf, emf[ , plotx2] >= x2lim[1] & emf[ , plotx2] <= x2lim[2])
    
    x1 <- emf[, plotx1]
    x2 <- emf[, plotx2]
    if (!is.numeric(x1))
        stop(paste("plotPlane: The variable", plotx1, "should be a numeric variable"))
    if (!is.numeric(x2))
        stop(paste("plotPlane: The variable", plotx2, "should be a numeric variable"))
    x1range <- magRange(x1, 1.15)
    x2range <- magRange(x2, 1.15)
  
    if (missing(ylab)) ylab <- colnames(model$model)[1]
    if (missing(x1lab)) x1lab <- plotx1
    if (missing(x2lab)) x2lab <- plotx2

   
    if (missing(x2lim)) x2range <- magRange(x2, 1.15) else x2range <- magRange(x2lim, 1.15)

    ##TODO must double check effect of function predictors
    otherPredictors <- setdiff(varnames, c(plotx2, plotx1)) #remove x1 x2
    if (length(otherPredictors) > 0) {
        otherPredictorValues <- centralValues(emf[ , otherPredictors, drop = FALSE])
    }

    myPredict <- function(a,b){
        ndf <- data.frame(a, b) #ndf = new data frame
        colnames(ndf) <- c(plotx1, plotx2)
        if (length(otherPredictors) > 0) {
            ndf <- cbind(ndf, otherPredictorValues)
            colnames(ndf) <- c(plotx1, plotx2, otherPredictors)
        }
        if ("glm" %in% class(model)) {
            if ("mcreg" %in% class(model)) attr(ndf, "isCentered") <- TRUE
            predict(model, newdata = ndf, type = "response")
        } else {
            if ("mcreg" %in% class(model)) attr(ndf, "isCentered") <- TRUE
            predict(model, newdata = ndf)
        }
    }

    x1seq <- plotSeq(x1range, length.out = npp)
    x2seq <- plotSeq(x2range, length.out = npp)
    zplane <- outer(x1seq, x2seq, function(a, b) { myPredict(a,b) } )

    yrange <- magRange(c(zplane,y), 1.15)

    ## must block users from trying to send wrong messages to persp.
    ## this is brutish, but don't have time figure more eloquent way.
    dotargs <- list(...)
    dotargs[["xlim"]] <- NULL
    dotargs[["ylim"]] <- NULL
    dotargs[["zlim"]] <- NULL

    argList <- list(x1 = plotSeq(x1range, x1floor),
                      x2 = plotSeq(x2range, x2floor), y = yrange,
                      x1lab = x1lab, x2lab = x2lab, ylab = ylab)
    newArgList <- modifyList(dotargs, argList)
    
    res <- do.call("perspEmpty", newArgList)

    ##for arrows. NEEDS reworking to be more general
    if ("glm" %in% class(model)) {
        fits <-  predict(model, type = "response")
    }else{
        fits <-  fitted(model)
    }

    mypoints4 <- trans3d(x1, x2, fits, pmat = res)
    newy <- ifelse(fits < y, fits + 0.8 * (y - fits),
                   fits + 0.8 * (y - fits))
    mypoints2s <- trans3d(x1, x2, newy, pmat = res)
    if (drawArrows){
        lmp <- length(fits)
        if (length(acol) > 1 && length(acol) < lmp) acol <- rep(acol, length.out = lmp)
        if (length(alty) > 1 && length(alty) < lmp) acol <- rep(acol, length.out = lmp)
        if (length(alwd) > 1 && length(alwd) < lmp) alwd <- rep(alwd, length.out = lmp)
        if (length(alength) > 1 && length(alength) < lmp) alwd <- rep(alength, length.out = lmp)
        arrows(mypoints4$x, mypoints4$y, mypoints2s$x, mypoints2s$y,
               col = acol, lty = alty, lwd = alwd, length = alength)
    }

    if (plotPoints){
        mypoints2 <- trans3d(x1, x2, y, pmat = res)
        lmp <- length(y)
        if (length(pch) > 1 && length(pch) < lmp) pch <- rep(pch, length.out = lmp)
        if (length(pcol) > 1 && length(pcol) < lmp) pcol <- rep(pcol, length.out = lmp)
        if (length(plwd) > 1 && length(plwd) < lmp) plwd <- rep(plwd, length.out = lmp)
        if (length(pcex) > 1 && length(pcex) < lmp) pcex <- rep(pcex, length.out = lmp)

        points(mypoints2, pch = pch, col = pcol, lwd = plwd, cex = pcex)
    }

    ##recycles parameters as necessary
    if (length(llwd) > 1 && length(llwd) < npp) llwd <- rep(llwd, length.out = npp)
    if (length(lcol) > 1 && length(lcol) < npp) lcol <- rep(lcol, length.out = npp)
    if (length(llty) > 1 && length(llty) < npp) llty <- rep(llty, length.out = npp)

    for (i in 1:npp) {
        lines(trans3d(x1seq[i], x2seq, zplane[i, ], pmat = res),
              lwd = llwd, col = lcol, lty= llty)
    }

    for (j in 1:npp) {
        lines(trans3d(x1seq, x2seq[j], zplane[, j], pmat = res),
              lwd = llwd, col = lcol, lty = llty)
    }


    if (!missing(linesFrom)) {
        dataSplits <- split(linesFrom$newdata, f = linesFrom$newdata[[linesFrom$call[["modx"]]]])
        drawLine <- function(nd, mycol, mylty){
            lines(trans3d(nd[[plotx1]], nd[[plotx2]], nd$fit, pmat=res), col = mycol, lwd = lflwd, lty = mylty)
        }
        mapply(drawLine, dataSplits, linesFrom$col, linesFrom$lty)
    }

    retval <- list(res = res, call = cl, "x1seq" = x1seq, "x2seq" = x2seq, "zplane" = zplane)
    class(retval) <- "rockchalk3d"
    invisible(retval)
}
NULL


## New for version 1.7.x

##' Superimpose regression lines on a plotted plane
##'
##' The examples will demonstrate the intended usage.
##'
##' From an educational stand point, the objective is to assist with the
##' student's conceptualization of the two and three dimensional regression
##' relationships.
##' @param to a 3d plot object produced by plotPlane
##' @param from output from a plotSlopes or plotCurves function (class="rockchalk")
##' @param col color of plotted lines (default: "red")
##' @param lwd line width of added lines (default: 2)
##' @param lty line type of added lines (default: 1)
##' @return NULL, nothing, nicht, nada.
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @export
##' @example inst/examples/addLines-ex.R

addLines <- function(to = NULL, from = NULL, col, lwd = 2, lty = 1){
    if (!inherits(from, "rockchalk"))
        stop("addLines: from must be an output object from plotSlopes or plotCurves, of class rockchalk")
    if (!inherits(to, "rockchalk3d"))
        stop("addLines: to must be a 3d plot object created by plotPlane or such")
    if ( !(from$call[["modx"]] %in% to$call[["plotx1"]] | from$call[["modx"]] %in% to$call[["plotx2"]]) )
        stop("Mismatched plotPlanes and plotSlopes objects")

    dataSplits <- split(from$newdata, f = from$newdata[[from$call[["modx"]]]])
    if(missing(col)) {
        col = from$col
    } else {
        if (length(col) < length(dataSplits)) col <- rep(col, length.out = length(dataSplits))
    }

    if (length(lwd) < length(dataSplits)) lwd <- rep(lwd, length.out = length(dataSplits))
    i <- 0
    for (i in seq_along(dataSplits)) {
        nd <- dataSplits[[i]]
        lines(trans3d(nd[[to$call$plotx1]], nd[[to$call$plotx2]],
                      nd$fit, pmat = to$res), col = col[i], lwd = lwd[i], lty = lty);
    }
    NULL
}



