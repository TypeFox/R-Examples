### File R/sisalTable.R
### This file is part of the sisal package for R.
###
### Copyright (C) 2015 Aalto University
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### A copy of the GNU General Public License is available at
### http://www.r-project.org/Licenses/

## Returns a grey value (as an R color string) that has largest
## minimum distance (best contrast) to bgLumas (numeric where 0 is
## black and 1 is white).
fgGrey <- function(bgLumas) {
    distWhite <- min(1 - bgLumas)
    distBlack <- min(bgLumas)
    meanLuma <- mean(bgLumas)
    distMean <- min(abs(meanLuma - bgLumas))
    switch(which.max(c(distBlack, distWhite, distMean)),
           "black", "white", rgb(meanLuma, meanLuma, meanLuma))
}

## A "standard" way of creating two-part names for viewports / grobs
twoPartName <- function(part1, part2) {
    paste(format(part1, scientific = FALSE, trim = TRUE, justify = "none"),
          format(part2, scientific = FALSE, trim = TRUE, justify = "none"),
          sep=".")
}
## when both part1 and part2 are already character vectors
twoPartNameChar <- function(part1, part2) {
    paste(part1, part2, sep=".")
}

### Drawing time auxiliary functions for sisalTable
leftRight.sisalTable <- function(x) {
    X <- convertUnit(x[["x"]], unitTo = "inches", axisFrom = "x",
                     typeFrom = "location", valueOnly = TRUE)
    hjust <- x[["internals"]][["hjust"]]
    width <- convertUnit(x[["width"]], unitTo = "inches", axisFrom = "x",
                         typeFrom = "dimension", valueOnly = TRUE)
    left <- X - width * hjust
    right <- left + width
    c(left, right)
}
bottomTop.sisalTable <- function(x) {
    Y <- convertUnit(x[["y"]], unitTo = "inches", axisFrom = "y",
                     typeFrom = "location", valueOnly = TRUE)
    vjust <- x[["internals"]][["vjust"]]
    height <- convertUnit(x[["height"]], unitTo = "inches", axisFrom = "y",
                          typeFrom = "dimension", valueOnly = TRUE)
    bottom <- Y - height * vjust
    top <- bottom + height
    c(bottom, top)
}

### S3 methods for sisalTable

## According to "frame" vignette, widthDetails and heightDetails
## should use absolute.size() which returns "null" units if the grob
## does not have an absolute size. I tried absolute.size(x$width),
## absolute.size(x$height), but that only gave an obscure error message:
##
## Error in UseMethod("absolute.units") :
##   no applicable method for 'absolute.units' applied to an object of class "unit"
## > traceback()
## 19: FUN(X[[1L]], ...)
## 18: lapply(unit, absolute.units)
## 17: absolute.units.unit.list(unit$arg1)
## 16: absolute.units(unit$arg1)
## 15: unit.arithmetic("sum", absolute.units(unit$arg1))
## 14: absolute.units.unit.arithmetic(unit)
## 13: absolute.units(unit)
## 12: absolute.size(x$width)
## 11: widthDetails.sisalTable(x)
## 10: widthDetails(x)
##
## The error message does not make sense. The (partial) traceback
## isn't too helpful either. So, let's use convertUnit(unitTo =
## "inches") instead.
##
## * Advantage: gives the size of the grob, even if the size is not
##   absolute (fixed). Seems reasonable to me.
##
## * Downside: does not work as discussed in the vignette. The
##   vignette may be outdated, though. For example, contrary to what
##   is stated in the vignette, neither widthDetails.rect() nor
##   heightDetails.rect() uses absolute.size(), which gets very little
##   use anywhere. Makes me wonder if there is a problem in
##   absolute.size(), or if I have overlooked something.
##
## This was checked on R trunk SVN revision 61876.
widthDetails.sisalTable <- function(x, ...) {
    convertUnit(x[["width"]], unitTo = "inches", axisFrom = "x",
                typeFrom = "dimension")
}
heightDetails.sisalTable <- function(x, ...) {
    convertUnit(x[["height"]], unitTo = "inches", axisFrom = "y",
                typeFrom = "dimension")
}
## xDetails and yDetails methods of sisalTable are vectorized
## w.r.t. theta, but it may not be possible to take advantage of this
## through grid
xDetails.sisalTable <- function(x, theta, ...) {
    leftRight <- leftRight.sisalTable(x)
    bottomTop <- bottomTop.sisalTable(x)
    width <- leftRight[2] - leftRight[1]
    height <- bottomTop[2] - bottomTop[1]
    X <- (leftRight[1] + leftRight[2]) / 2
    Y <- (bottomTop[1] + bottomTop[2]) / 2
    rGrob <- rectGrob(x = X, y = Y, width = width, height = height,
                      default.units = "inches")
    do.call("unit.c", lapply(theta, function (x) xDetails(rGrob, x)))
}
yDetails.sisalTable <- function(x, theta, ...) {
    leftRight <- leftRight.sisalTable(x)
    bottomTop <- bottomTop.sisalTable(x)
    width <- leftRight[2] - leftRight[1]
    height <- bottomTop[2] - bottomTop[1]
    X <- (leftRight[1] + leftRight[2]) / 2
    Y <- (bottomTop[1] + bottomTop[2]) / 2
    rGrob <- rectGrob(x = X, y = Y, width = width, height = height,
                      default.units = "inches")
    do.call("unit.c", lapply(theta, function (x) yDetails(rGrob, x)))
}
makeContent.sisalTable <- function(x, ...) {
    x2 <- x
    ## Edit content labels
    internals <- x[["internals"]]
    cellSide <- convertWidth(internals[["cellSide"]],
                             unitTo = "inches", valueOnly = TRUE)
    nRows <- internals[["nRows"]]
    nCols <- internals[["nCols"]]
    nResize <- internals[["nResizeLabels"]]
    if (nResize > 0) {
        labelFS <- cellSide / internals[["labelDim"]]
        labelsNotNA <- internals[["labelsNotNA"]]
        nLabels <- length(labelsNotNA)
        textNames <- internals[["textNames"]]
        labelPars <- internals[["labelPars"]]
        for (i in seq_along(labelPars)) {
            labelPars[[i]][["fontsize"]] <- labelFS
        }
        labelMeasurements <-
            dynTextGrob(internals[["labels"]], adjustJust = TRUE,
                        gp = labelPars, resize = FALSE,
                        takeMeasurements = TRUE)
        labelJust <- numeric(nLabels)
        justTmp <- labelMeasurements[["vjust"]]
        labelJust[labelsNotNA] <- justTmp[1:nResize]
        linearRowCol <- 0
        for (Col in seq_len(nCols)) {
            ## Actual contents of the table
            for (Row in seq_len(nRows)) {
                labelsMod <- linearRowCol %% nLabels + 1
                if (labelsNotNA[labelsMod]) {
                    thisName <- textNames[linearRowCol + 1]
                    labelPar <- gparToList(getGrob(x, thisName)[["gp"]])
                    labelPar[["fontsize"]] <- labelFS
                    x2 <- editGrob(x2, thisName,
                                   gp = do.call("gpar", labelPar),
                                   vjust = labelJust[labelsMod])
                }
                linearRowCol <- linearRowCol + 1
            }
        }
    }
    if (internals[["resizeAxes"]]) {
        ## Edit x axis labels
        axesFS <- NULL
        xAxisSide <- internals[["xAxisSide"]]
        if (!is.null(xAxisSide)) {
            xAxisSide <- convertHeight(xAxisSide,
                                       unitTo = "inches", valueOnly = TRUE)
            axesFS <- min(cellSide / internals[["xAxisWidth"]],
                          xAxisSide / internals[["xAxisHeight"]])
            xAxisNotNA <- internals[["xAxisNotNA"]]
            axesPar <- internals[["axesPar"]]
            axesPar[["fontsize"]] <- axesFS
            axesPar <- do.call("gpar", axesPar)
            xNames <- internals[["xNames"]]
            xMeasurements <-
                dynTextGrob(internals[["xAxisLabels"]], resize = FALSE,
                            gp = axesPar, rotJust = TRUE,
                            rot = internals[["xAxisRot"]],
                            takeMeasurements = TRUE,
                            just = internals[["xAxisJust"]], adjustJust = TRUE)
            hjust <- numeric(nCols)
            vjust <- numeric(nCols)
            hjust[xAxisNotNA] <- xMeasurements[["hjust"]]
            vjust[xAxisNotNA] <- xMeasurements[["vjust"]]
            xShift <- numeric(nCols)
            yShift <- numeric(nCols)
            xShift[xAxisNotNA] <- xMeasurements[["xShift"]]
            yShift[xAxisNotNA] <- xMeasurements[["yShift"]]
            nXAxis <- length(xAxisNotNA)
            for (Col in seq_len(nCols)[xAxisNotNA]) {
                xMod <- (Col - 1) %% nXAxis + 1
                thisName <- xNames[Col]
                thisGrob <- getGrob(x, thisName)
                x2 <- editGrob(x2, thisName, gp = axesPar,
                               x=thisGrob[["x"]] + unit(xShift[xMod],"inches"),
                               y=thisGrob[["y"]] + unit(yShift[xMod],"inches"),
                               hjust = hjust[xMod], vjust = vjust[xMod])
            }
        }
        ## Edit y axis labels
        yAxisSide <- internals[["yAxisSide"]]
        if (!is.null(yAxisSide)) {
            if (is.null(axesFS)) {
                yAxisSide <- convertWidth(yAxisSide,
                                          unitTo = "inches", valueOnly = TRUE)
                axesFS <- min(cellSide / internals[["yAxisHeight"]],
                              yAxisSide / internals[["yAxisWidth"]])
                axesPar <- internals[["axesPar"]]
                axesPar[["fontsize"]] <- axesFS
                axesPar <- do.call("gpar", axesPar)
            }
            yAxisNotNA <- internals[["yAxisNotNA"]]
            yNames <- internals[["yNames"]]
            yMeasurements <-
                dynTextGrob(internals[["yAxisLabels"]], resize = FALSE,
                            gp = axesPar, rotJust = TRUE,
                            rot = internals[["yAxisRot"]],
                            takeMeasurements = TRUE,
                            just = internals[["yAxisJust"]], adjustJust = TRUE)
            hjust <- numeric(nRows)
            vjust <- numeric(nRows)
            hjust[yAxisNotNA] <- yMeasurements[["hjust"]]
            vjust[yAxisNotNA] <- yMeasurements[["vjust"]]
            xShift <- numeric(nRows)
            yShift <- numeric(nRows)
            xShift[yAxisNotNA] <- yMeasurements[["xShift"]]
            yShift[yAxisNotNA] <- yMeasurements[["yShift"]]
            nYAxis <- length(yAxisNotNA)
            for (Row in seq_len(nRows)[yAxisNotNA]) {
                yMod <- (Row - 1) %% nYAxis + 1
                thisName <- yNames[Row]
                thisGrob <- getGrob(x, thisName)
                x2 <- editGrob(x2, thisName, gp = axesPar,
                               x=thisGrob[["x"]] + unit(xShift[yMod],"inches"),
                               y=thisGrob[["y"]] + unit(yShift[yMod],"inches"),
                               hjust = hjust[yMod], vjust = vjust[yMod])
            }
        }
    }
    ## Final edited grobTree
    x2
}
validDetails.sisalTable <- function(x, ...) {
    if (!inherits(x, "gTree")) {
        stop("a sisalTable object must be a gTree")
    }
    slotNames <- names(x)
    requiredNames <- c("width", "height", "vpLayout", "x", "y", "just",
                       "clip", "internals", "name", "gp", "vp")

    ## Testing that all slots are present
    if (!all(requiredNames %in% slotNames)) {
        stop("slots missing")
    }

    ## Testing 'x', 'y', 'width', 'height', 'clip', 'just'
    X <- x[["x"]]
    Y <- x[["y"]]
    width <- x[["width"]]
    height <- x[["height"]]
    if (!inherits(X, "unit") || !inherits(Y, "unit") ||
        !inherits(width, "unit") || !inherits(height, "unit")) {
        stop("'x', 'y', 'width' and 'height' must be grid units")
    }
    if (length(X) != 1 || length(Y) != 1 ||
        length(width) != 1 || length(height) != 1) {
        stop("'x', 'y', 'width' and 'height' must have one element each")
    }
    clip <- x[["clip"]]
    if (!(length(clip) == 1 &&
          ((is.character(clip) && clip %in% c("on", "inherit", "off")) ||
           (is.logical(clip) && is.finite(clip))))) {
        stop("'clip' must be TRUE, FALSE, \"on\", \"inherit\" or \"off\"")
    }
    just <- x[["just"]]
    if (!(((is.numeric(just) && all(is.finite(just))) ||
           is.character(just)) && length(just) > 0 && length(just) < 3)) {
        stop("'just' must be a finite numeric or character vector of length 1 or 2")
    }

    ## Testing internal variables
    internals <- x[["internals"]]
    if (!is.list(internals)) {
        stop(gettextf("'%s' must be a list", 'x[["internals"]]',
                      domain = "R-sisal"), domain = NA)
    }
    varNames <-
        c("hjust", "vjust", "nResizeLabels", "resizeAxes", "labelDim",
          "xAxisWidth", "xAxisHeight", "yAxisWidth", "yAxisHeight", "cellSide",
          "xAxisSide", "yAxisSide", "nRows", "nCols", "labelsNotNA",
          "xAxisNotNA", "yAxisNotNA", "labels", "labelPars", "axesPar",
          "xAxisLabels", "yAxisLabels", "xAxisRot", "yAxisRot", "xAxisJust",
          "yAxisJust", "textNames", "xNames", "yNames", "childVps")
    if (!all(varNames %in% names(internals))) {
        stop(gettextf("'%s' must have the required variables",
                      'x[["internals"]]', domain = "R-sisal"), domain = NA)
    }
    hjust <- internals[["hjust"]]
    vjust <- internals[["vjust"]]
    nResize <- internals[["nResizeLabels"]]
    resizeAxes <- internals[["resizeAxes"]]
    labelDim <- internals[["labelDim"]]
    xAxisWidth <- internals[["xAxisWidth"]]
    xAxisHeight <- internals[["xAxisHeight"]]
    yAxisWidth <- internals[["yAxisWidth"]]
    yAxisHeight <- internals[["yAxisHeight"]]
    cellSide <- internals[["cellSide"]]
    xAxisSide <- internals[["xAxisSide"]]
    yAxisSide <- internals[["yAxisSide"]]
    nRows <- internals[["nRows"]]
    nCols <- internals[["nCols"]]
    labelsNotNA <- internals[["labelsNotNA"]]
    labels <- internals[["labels"]]
    labelPars <- internals[["labelPars"]]
    xAxisNotNA <- internals[["xAxisNotNA"]]
    yAxisNotNA <- internals[["yAxisNotNA"]]
    xAxisLabels <- internals[["xAxisLabels"]]
    yAxisLabels <- internals[["yAxisLabels"]]
    xAxisRot <- internals[["xAxisRot"]]
    yAxisRot <- internals[["yAxisRot"]]
    xAxisJust <- internals[["xAxisJust"]]
    yAxisJust <- internals[["yAxisJust"]]
    textNames <- internals[["textNames"]]
    xNames <- internals[["xNames"]]
    yNames <- internals[["yNames"]]
    stopifnot(is.numeric(hjust), is.numeric(vjust), is.numeric(nResize),
              is.logical(resizeAxes), is.numeric(labelDim),
              is.numeric(xAxisWidth), is.numeric(xAxisHeight),
              is.numeric(yAxisWidth), is.numeric(yAxisHeight),
              inherits(cellSide, "unit"), is.numeric(nRows), is.numeric(nCols),
              is.list(internals[["axesPar"]]), is.numeric(xAxisRot),
              is.numeric(yAxisRot),
              inherits(internals[["childVps"]], "vpList"))
    stopifnot(length(hjust) == 1, length(vjust) == 1, length(nResize) == 1,
              length(resizeAxes) == 1, length(labelDim) == 1,
              length(xAxisWidth) == 1, length(xAxisHeight) == 1,
              length(yAxisWidth) == 1, length(yAxisHeight) == 1,
              length(cellSide) == 1, length(nRows) == 1, length(nCols) == 1,
              length(xAxisRot) == 1, length(yAxisRot) == 1)
    stopifnot(is.finite(hjust), is.finite(vjust),
              is.finite(nResize) && nResize >= 0, is.finite(resizeAxes),
              is.finite(labelDim) && labelDim >= 0,
              is.finite(nRows) && nRows > 0, is.finite(nCols) && nCols > 0,
              is.finite(xAxisRot), is.finite(yAxisRot))
    stopifnot(!resizeAxes || is.null(xAxisSide) ||
              (is.finite(xAxisWidth) && xAxisWidth >= 0 &&
               is.finite(xAxisHeight) && xAxisHeight >= 0 &&
               inherits(xAxisSide, "unit") && length(xAxisSide) == 1 &&
               is.logical(xAxisNotNA) && length(xAxisNotNA) >= 1 &&
               is.finite(xAxisNotNA) && is.vector(xAxisLabels) &&
               length(xAxisLabels) >= 1 && is.character(xNames) &&
               length(xNames) == nCols))
    stopifnot(!resizeAxes || is.null(yAxisSide) ||
              (is.finite(yAxisWidth) && yAxisWidth >= 0 &&
               is.finite(yAxisHeight) && yAxisHeight >= 0 &&
               inherits(yAxisSide, "unit") && length(yAxisSide) == 1 &&
               is.logical(yAxisNotNA) && length(yAxisNotNA) >= 1 &&
               is.finite(yAxisNotNA) && is.vector(yAxisLabels) &&
               length(yAxisLabels) >= 1 && is.character(yNames) &&
               length(yNames) == nRows))
    stopifnot(length(xAxisJust) > 0,
              is.character(xAxisJust) ||
              (is.numeric(xAxisJust) && all(is.finite(xAxisJust))))
    stopifnot(length(yAxisJust) > 0,
              is.character(yAxisJust) ||
              (is.numeric(yAxisJust) && all(is.finite(yAxisJust))))
    if (nResize > 0) {
        stopifnot(is.logical(labelsNotNA), is.finite(labelsNotNA),
                  is.character(textNames), length(textNames) == nRows * nCols,
                  is.list(labelPars), length(labelPars) >= 2,
                  vapply(labelPars, is.list, logical(1)),
                  is.vector(labels), length(labels) >= 2)
    }
    x
}
editDetails.sisalTable <- function(x, specs) {
    x2 <- NextMethod("editDetails", x)
    slotNames <- names(specs)
    if ("internals" %in% slotNames) {
        stop("'internals' must not be edited")
    }
    if (any(c("width", "height", "x", "y",
              "just", "clip", "vpLayout") %in% slotNames)) {
        internals <- x2[["internals"]]
        parentVp <- viewport(x = x2[["x"]], y = x2[["y"]],
                             width = x2[["width"]], height = x2[["height"]],
                             just = x2[["just"]], clip = x2[["clip"]],
                             layout = x2[["vpLayout"]], name = "Parent")
        x2 <- editGrob(x2, childrenvp = vpTree(parent = parentVp,
                               children = internals[["childVps"]]))
        if ("just" %in% slotNames) {
            x2[["internals"]][c("hjust", "vjust")] <- justHV(x2[["just"]])
        }
    }
    x2
}

### Main functions

## Computes sizes of different graphical elements.  Returns those and
## suggested relative space allocations for the various parts of a
## table (based on computed sizes and target ratios specified as
## arguments).  All the gp* arguments should be (lists of) ordinary
## lists containing settings accepted by gpar().
tableMeasures <- function(xAxisLabels, yAxisLabels, labels, targetAxesRatio,
                          gpAxes = list(), gpLabels,
                          main, xlab, ylab, gpMain, gpLab, nRows, nCols,
                          targetMainRatio = 0.15, targetLabRatio = 0.15,
                          scaleLabels = FALSE, xlabExpand = 1, ylabExpand = 1,
                          mainExpand = 1, xAxisRot = 0, yAxisRot = 0,
                          xAxisJust, yAxisJust) {
    nX <- length(xAxisLabels)
    if (nX > 0) {
        tmp <- dynTextGrob(xAxisLabels, gp = do.call("gpar", gpAxes),
                           rot = xAxisRot, rotJust = TRUE,
                           just = xAxisJust, takeMeasurements = TRUE,
                           adjustJust = TRUE)
        xAxisWidth <- tmp[["sizingWidth"]]
        xAxisHeight <- tmp[["sizingHeight"]]
        xAxisNZ <- tmp[["nzDim"]]
    } else {
        xAxisWidth <- NA_real_
        xAxisHeight <- NA_real_
        xAxisNZ <- FALSE
    }
    nY <- length(yAxisLabels)
    if (nY > 0) {
        tmp <- dynTextGrob(yAxisLabels, gp = do.call("gpar", gpAxes),
                           rot = yAxisRot, rotJust = TRUE,
                           just = yAxisJust, takeMeasurements = TRUE,
                           adjustJust = TRUE)
        yAxisWidth <- tmp[["sizingWidth"]]
        yAxisHeight <- tmp[["sizingHeight"]]
        yAxisNZ <- tmp[["nzDim"]]
    } else {
        yAxisWidth <- NA_real_
        yAxisHeight <- NA_real_
        yAxisNZ <- FALSE
    }
    labelCex <- 0
    nLabels <- length(labels)
    for (i in seq_along(gpLabels)) {
        cexTemp <- gpLabels[[i]][["cex"]]
        if (is.null(cexTemp)) {
            labelCex <- max(1, labelCex)
        } else {
            labelCex <- max(cexTemp, labelCex)
        }
    }
    mfs <- dynTextGrob(c(labels, "x"), adjustJust = TRUE,
                       gp = c(gpLabels, list(list())), takeMeasurements = TRUE)
    labelDim <- max(unlist(mfs[c("sizingWidth", "sizingHeight")]))
    tmpSeq <- seq_len(nLabels)
    labelNZ <- mfs[["nzDim"]][tmpSeq]
    if (length(main) > 0) {
        tmp <- dynTextGrob(main, gp = gpMain, takeMeasurements = TRUE)
        maxMainWidth <- tmp[["sizingWidth"]]
        maxMainHeight <- tmp[["sizingHeight"]]
        mainNZ <- tmp[["nzDim"]]
    } else {
        maxMainWidth <- NA_real_
        maxMainHeight <- NA_real_
        mainNZ <- FALSE
    }
    if (length(xlab) > 0) {
        tmp <- dynTextGrob(xlab, gp = gpLab, takeMeasurements = TRUE)
        maxXlabWidth <- tmp[["sizingWidth"]]
        maxXlabHeight <- tmp[["sizingHeight"]]
        xlabNZ <- tmp[["nzDim"]]
    } else {
        maxXlabWidth <- NA_real_
        maxXlabHeight <- NA_real_
        xlabNZ <- FALSE
    }
    if (length(ylab) > 0) {
        tmp <- dynTextGrob(ylab, gp = gpLab, rot = 90, takeMeasurements = TRUE)
        maxYlabWidth <- tmp[["sizingWidth"]]
        maxYlabHeight <- tmp[["sizingHeight"]]
        ylabNZ <- tmp[["nzDim"]]
    } else {
        maxYlabWidth <- NA_real_
        maxYlabHeight <- NA_real_
        ylabNZ <- FALSE
    }
    targetAxesRatio2 <- targetAxesRatio * labelCex
    if (is.na(xAxisHeight)) {
        xAxisRatio <- Inf
    } else if (targetAxesRatio2 * xAxisWidth < labelDim) {
        xAxisRatio <- xAxisHeight * targetAxesRatio2 / labelDim
    } else {
        xAxisRatio <- xAxisHeight / xAxisWidth
    }
    if (is.na(yAxisWidth)) {
        yAxisRatio <- Inf
    } else if (targetAxesRatio2 * yAxisHeight < labelDim) {
        yAxisRatio <- yAxisWidth * targetAxesRatio2 / labelDim
    } else {
        yAxisRatio <- yAxisWidth / yAxisHeight
    }
    xFs <- min(1 / xAxisWidth, xAxisRatio / xAxisHeight)
    yFs <- min(1 / yAxisHeight, yAxisRatio / yAxisWidth)
    xAxisRatio <- min(xAxisRatio, xAxisRatio * yFs / xFs, na.rm = TRUE)
    yAxisRatio <- min(yAxisRatio, yAxisRatio * xFs / yFs, na.rm = TRUE)
    if (!is.na(xFs) || !is.na(yFs)) {
        xyFs <- min(xFs, yFs, na.rm = TRUE)
        if (scaleLabels && !is.na(labelDim) &&
            xyFs < targetAxesRatio2 / labelDim) {
            fontScaleIn <- xyFs / (targetAxesRatio2 / labelDim)
            labelDim <- labelDim / fontScaleIn
        }
    } else {
        xyFs <- NA_real_
    }

    axesCex <- gpAxes[["cex"]]
    if (is.null(axesCex)) {
        axesCex <- 1
    }
    mainRatio <- min(min(1, mainExpand) * nCols * maxMainHeight / maxMainWidth,
                     mainExpand * max(targetMainRatio * nRows,
                                      maxMainHeight * labelCex / labelDim,
                                      maxMainHeight * axesCex * xyFs,
                                      na.rm = TRUE), na.rm = TRUE)
    xlabRatio <- min(min(1, xlabExpand) * nCols * maxXlabHeight / maxXlabWidth,
                     xlabExpand * max(targetLabRatio * min(nRows, nCols),
                                      maxXlabHeight * labelCex / labelDim,
                                      maxXlabHeight * axesCex * xyFs,
                                      na.rm = TRUE), na.rm = TRUE)
    ylabRatio <- min(min(1, ylabExpand) * nRows * maxYlabWidth / maxYlabHeight,
                     ylabExpand * max(targetLabRatio * min(nRows, nCols),
                                      maxYlabWidth * labelCex / labelDim,
                                      maxYlabWidth * axesCex * xyFs,
                                      na.rm = TRUE), na.rm = TRUE)
    ylabScale <- min(nRows / maxYlabHeight, ylabRatio / maxYlabWidth) /
        ylabExpand
    xlabScale <- min(nCols / maxXlabWidth, xlabRatio / maxXlabHeight) /
        xlabExpand
    if (!is.na(ylabScale) && !is.na(xlabScale)) {
        if (ylabScale > xlabScale) {
            ylabRatio <- xlabScale / ylabScale * ylabRatio
        } else {
            xlabRatio <- ylabScale / xlabScale * xlabRatio
        }
    }

    ## * *Ratio: Ratio of space reserved (along a dimension) for
    ## various graphical objects to length of the sides of each
    ## (square) cell in the table

    ## * *Dim, *Width, *Height: Maximum normalized size (Dim: max of
    ## width and height) of a graphical object. Normalized means size
    ## of the graphical object in points divided by the specified font
    ## size.
    list(xAxisRatio = xAxisRatio,        # x axis labels (one for each column)
         yAxisRatio = yAxisRatio,        # y axis labels (one for each row)
         labelDim = labelDim,            # labels in the table
         labelNZ = labelNZ,
         xAxisWidth = xAxisWidth,        # x axis labels
         xAxisHeight = xAxisHeight,      # x axis labels
         xAxisNZ = xAxisNZ,
         yAxisWidth = yAxisWidth,        # y axis labels
         yAxisHeight = yAxisHeight,      # y axis labels
         yAxisNZ = yAxisNZ,
         mainRatio = mainRatio,          # main title
         mainWidth = maxMainWidth,       # main title
         mainHeight = maxMainHeight,     # main title
         mainNZ = mainNZ,
         xlabRatio = xlabRatio,          # x label (one across all cols)
         ylabRatio = ylabRatio,          # y label (one across all rows)
         xlabWidth = maxXlabWidth,       # x label
         xlabHeight = maxXlabHeight,     # x label
         ylabWidth = maxYlabWidth,       # y label
         ylabHeight = maxYlabHeight,     # y label
         xlabNZ = xlabNZ,
         ylabNZ = ylabNZ)
}

## Convert unit() object x to inches (value only) in the current
## context.  Uses both axisFrom = "x" and axisFrom = "y" and takes the
## minimum.
minConvertXY <- function(x) {
    min(convertUnit(x, unitTo = "inches",
                    axisFrom = "x", typeFrom = "dimension",
                    valueOnly = TRUE),
        convertUnit(x, unitTo = "inches",
                    axisFrom = "y", typeFrom = "dimension",
                    valueOnly = TRUE))
}

sisalTable <- function(labels = matrix(seq_len(12), 3, 4),
                       nRows = NROW(labels), nCols = NCOL(labels),
                       bg = sample(colors(), nRows * nCols, replace=TRUE),
                       stripeCol = NULL,
                       fg = NULL, naFill = "white", naStripes = "grey50",
                       main = NULL, xlab = NULL, ylab = NULL,
                       xAxisLabels = NULL, yAxisLabels = NULL,
                       draw = TRUE, outerRect = TRUE, innerLines = TRUE,
                       nStripes = 7, stripeRot = 45, stripeWidth = 0.2,
                       stripeScale = 0.95, resizeText = TRUE,
                       resizeTable = TRUE, resizeMain = resizeText,
                       resizeLab = resizeText, resizeAxes = resizeText,
                       resizeLabels = resizeTable && resizeText,
                       x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                       width = unit(0.97, "npc"), height = unit(0.97, "npc"),
                       default.units = "npc", just = "center",
                       clip = "inherit", xAxisRot = 0, yAxisRot = 0,
                       xAxisJust = c(0.5, 1), xAxisX = 0.5, xAxisY = 1,
                       yAxisJust = c(1, 0.5), yAxisX = 1, yAxisY = 0.5,
                       mainMargin = if(resizeMain) 0.15 else unit(8, "points"),
                       xlabMargin = if (resizeLab) 0.1 else unit(5, "points"),
                       ylabMargin = if (resizeLab) 0.1 else unit(5, "points"),
                       axesMargin = if (resizeAxes) 0.1 else unit(5, "points"),
                       axesSize = 0.8, forceAxesSize = FALSE, mainSize = 1,
                       xlabSize = 1, ylabSize = 1,
                       mainPar = gpar(fontface="bold", fontsize=14),
                       labPar = gpar(fontface="plain", fontsize=14),
                       labelPars = gpar(fontsize = 20, cex = 0.6),
                       axesPar = gpar(fontsize = 10),
                       rectPar = gpar(), linePar = gpar(),
                       name = NULL, gp = NULL, vp = NULL) {

    ## Argument checking
    stopifnot(is.numeric(xAxisX), length(xAxisX) == 1,
              is.finite(xAxisX), is.numeric(xAxisY), length(xAxisY) == 1,
              is.finite(xAxisY), is.numeric(yAxisX), length(yAxisX) == 1,
              is.finite(yAxisX), is.numeric(yAxisY), length(yAxisY) == 1,
              is.finite(yAxisY))
    stopifnot(identical(resizeMain, TRUE) || identical(resizeMain, FALSE),
              identical(resizeLab, TRUE) || identical(resizeLab, FALSE),
              identical(resizeLabels, TRUE) || identical(resizeLabels, FALSE),
              identical(resizeAxes, TRUE) || identical(resizeAxes, FALSE),
              identical(resizeTable, TRUE) || identical(resizeTable, FALSE),
              identical(draw, TRUE) || identical(draw, FALSE),
              identical(outerRect, TRUE) || identical(outerRect, FALSE),
              identical(innerLines, TRUE) || identical(innerLines, FALSE),
              identical(forceAxesSize, TRUE) ||
              identical(forceAxesSize, FALSE))
    stopifnot(length(just) > 0,
              is.character(just) || (is.numeric(just) && all(is.finite(just))))
    stopifnot(length(xAxisJust) > 0,
              is.character(xAxisJust) ||
              (is.numeric(xAxisJust) && all(is.finite(xAxisJust))))
    stopifnot(length(yAxisJust) > 0,
              is.character(yAxisJust) ||
              (is.numeric(yAxisJust) && all(is.finite(yAxisJust))))
    stopifnot(length(just) > 0,
              is.character(just) || (is.numeric(just) && all(is.finite(just))))
    stopifnot(is.null(name) ||
              (is.character(name) && length(name) == 1 && !is.na(name)))
    stopifnot(inherits(mainPar, "gpar"), inherits(labPar, "gpar"),
              inherits(labelPars, "gpar") ||
              (is.list(labelPars) &&
               all(vapply(labelPars, inherits, TRUE, what = "gpar"))),
              inherits(axesPar, "gpar"), inherits(rectPar, "gpar"),
              inherits(linePar, "gpar"))
    stopifnot(is.numeric(axesSize), length(axesSize) == 1, is.finite(axesSize),
              axesSize > 0, is.numeric(mainSize), length(mainSize) == 1,
              is.finite(mainSize), mainSize > 0, is.numeric(xlabSize),
              length(xlabSize) == 1, is.finite(xlabSize), xlabSize > 0,
              is.numeric(ylabSize), length(ylabSize) == 1, is.finite(ylabSize),
              ylabSize > 0)
    stopifnot(is.numeric(nRows), length(nRows) == 1,
              is.finite(nRows), nRows > 0, round(nRows) == nRows,
              is.numeric(nCols), length(nCols) == 1, is.finite(nCols),
              nCols > 0, round(nCols) == nCols)
    stopifnot(is.numeric(xAxisRot), length(xAxisRot) == 1, is.finite(xAxisRot),
              is.numeric(yAxisRot), length(yAxisRot) == 1, is.finite(yAxisRot))
    stopifnot(is.numeric(stripeRot), length(stripeRot) == 1,
              is.finite(stripeRot), is.numeric(nStripes),
              length(nStripes) == 1,
              is.finite(nStripes), nStripes >= 1, round(nStripes) == nStripes,
              is.numeric(stripeWidth), length(stripeWidth) == 1,
              is.finite(stripeWidth), stripeWidth > 0, stripeWidth < 1,
              is.numeric(stripeScale), length(stripeScale) == 1,
              is.finite(stripeScale), stripeScale >= 0, stripeScale <= 1)

    allNA <- function (x) {
        all(suppressWarnings(is.na(x))) || identical(x, "")
    }
    linearizeDf <- function(x) {
        do.call("c",
                lapply(x,
                       function(x) {
                           if (is.factor(x)) {
                               as.list(levels(x)[x])
                           } else {
                               as.list(x)
                           }
                       }))
    }

    ## resizeLabels doesn't count in anyResize
    anyResize <- any(resizeTable, resizeMain, resizeLab, resizeAxes)
    if (!draw && dev.cur() == 1) {
        ## Open dummy device => avoid opening a window when using
        ## convertUnit(). As far as I can see, the unit conversions
        ## used here are independent of device parameters and thus
        ## should not be confused by the dummy device.
        tFile <- tempfile()
        pdf(file = tFile) # pdf device is always available
        tDev <- dev.cur()
        on.exit(dev.off(tDev))
        on.exit(unlink(tFile), add = TRUE)
    }

    nCells <- nRows * nCols
    if (is.data.frame(bg)) {
        bg2 <- linearizeDf(bg)
    } else {
        bg2 <- bg
    }
    if (length(bg2) == 0) {
        bg2 <- rep.int(NA_character_, nCells)
    } else {
        bg2 <- rep_len(bg2, nCells)
    }
    if (is.data.frame(stripeCol)) {
        stripeCol2 <- linearizeDf(stripeCol)
    } else {
        stripeCol2 <- stripeCol
    }
    if (length(stripeCol2) == 0) {
        stripeCol2 <- rep.int(NA_character_, nCells)
    } else {
        stripeCol2 <- rep_len(stripeCol2, nCells)
    }

    ## Conversion coefficients from RGB to brightness / luminance / luma...
    ## Sources referenced on 2015-10-09.
    ## http://stackoverflow.com/q/596216
    ## https://en.wikipedia.org/wiki/Luma_%28video%29
    LUMA.CONVERSION <- c(0.2126, 0.7152, 0.0722) # R, G, B

    ## Get RGB values of bg colors (checks that bg colors are valid)
    bgNA <- is.na(bg2)
    bgRGB <- t(col2rgb(bg2[!bgNA]))
    stripeNA <- is.na(stripeCol2)
    bgLuma <- numeric(nCells)
    bgAndStripes <- which(!bgNA & !stripeNA)
    onlyBg <- !bgNA & stripeNA
    foo <- stripeNA[!bgNA]
    bgRGBBoth <- bgRGB[!foo, , drop=FALSE]
    bgRGBOnly <- bgRGB[foo, , drop=FALSE]
    stripeRGB <- t(col2rgb(stripeCol2[bgAndStripes]))
    luma1 <- drop(stripeRGB %*% LUMA.CONVERSION) / 255
    dist1 <- abs(luma1 - 0.5)
    luma2 <- drop(bgRGBBoth %*% LUMA.CONVERSION) / 255
    dist2 <- abs(luma2 - 0.5)
    dist1Larger <- dist1 > dist2
    dist2Larger <- !dist1Larger
    bgLuma[bgAndStripes[dist1Larger]] <- luma1[dist1Larger]
    bgLuma[bgAndStripes[dist2Larger]] <- luma2[dist2Larger]
    bgLuma[onlyBg] <- drop(bgRGBOnly %*% LUMA.CONVERSION / 255)
    bgLuma[bgNA] <-  drop(c(stripeWidth, 1 - stripeWidth) %*%
                          (t(col2rgb(c(naStripes, naFill))) %*%
                           LUMA.CONVERSION)) / 255
    if (is.data.frame(fg)) {
        fg2 <- linearizeDf(fg)
    } else {
        fg2 <- fg
    }
    if (length(fg2) == 0) {
        fg2 <- rep.int("black", nCells)
        fg2[bgLuma < 0.5] <- "white"
    } else if (length(fg2) > nCells) {
        fg2 <- fg2[seq_len(nCells)]
    }
    nFg <- length(fg2)
    if (is.data.frame(labels)) {
        labels2 <- linearizeDf(labels)
    } else {
        labels2 <- labels
    }
    if (length(labels2) == 0) {
        labelsNA <- TRUE
        nLabels <- 1
        labels2 <- character(0)
    } else if (length(labels2) > nCells) {
        labels2 <- labels2[seq_len(nCells)]
    }
    ## labels2,xAxisLabels2, yAxisLabels2 will be indexed with `[[`
    if (is.expression(labels2)) {
        labelTmp <- vector(mode = "list", length = length(labels2))
        for (i in seq_along(labelTmp)) {
            labelTmp[[i]] <- labels2[i]
        }
        labels2 <- labelTmp
    } else if (!is.vector(labels2)) {
        labels2 <- as.list(labels2)
    }
    if (length(labels2) > 0) {
        labelsNA <- vapply(labels2, allNA, TRUE, USE.NAMES = FALSE)
        nLabels <- length(labels2)
    }
    if (length(xAxisLabels) == 0) {
        xAxisNA <- TRUE
        nXAxis <- 1
        have.xAxisLabels <- FALSE
        xAxisLabels2 <- character(0)
    } else if (length(xAxisLabels) > nCols) {
        xAxisLabels2 <- xAxisLabels[seq_len(nCols)]
    } else {
        xAxisLabels2 <- xAxisLabels
    }
    if (is.expression(xAxisLabels2)) {
        labelTmp <- vector(mode = "list", length = length(xAxisLabels2))
        for (i in seq_along(labelTmp)) {
            labelTmp[[i]] <- xAxisLabels2[i]
        }
        xAxisLabels2 <- labelTmp
    } else if (!is.vector(xAxisLabels2)) {
        xAxisLabels2 <- as.list(xAxisLabels2)
    }
    if (length(xAxisLabels2) > 0) {
        xAxisNA <- vapply(xAxisLabels2, allNA, TRUE, USE.NAMES = FALSE)
        nXAxis <- length(xAxisLabels2)
        have.xAxisLabels <- !all(xAxisNA)
    }
    if (length(yAxisLabels) == 0) {
        yAxisNA <- TRUE
        nYAxis <- 1
        have.yAxisLabels <- FALSE
        yAxisLabels2 <- character(0)
    } else if (length(yAxisLabels) > nRows) {
        yAxisLabels2 <- yAxisLabels[seq_len(nRows)]
    } else {
        yAxisLabels2 <- yAxisLabels
    }
    if (is.expression(yAxisLabels2)) {
        labelTmp <- vector(mode = "list", length = length(yAxisLabels2))
        for (i in seq_along(labelTmp)) {
            labelTmp[[i]] <- yAxisLabels2[i]
        }
        yAxisLabels2 <- labelTmp
    } else if (!is.vector(yAxisLabels2)) {
        yAxisLabels2 <- as.list(yAxisLabels2)
    }
    if (length(yAxisLabels2) > 0) {
        yAxisNA <- vapply(yAxisLabels2, allNA, TRUE, USE.NAMES = FALSE)
        nYAxis <- length(yAxisLabels2)
        have.yAxisLabels <- !all(yAxisNA)
    }
    if (inherits(labelPars, "gpar")) {
        labelPars2 <- list(gparToList(labelPars))
        nLabelPars <- 1
    } else if (length(labelPars) == 0) {
        labelPars2 <- list(list())
        nLabelPars <- 1
    } else if (length(labelPars) > nCells) {
        labelPars2 <- lapply(labelPars[seq_len(nCells)], "gparToList")
        nLabelPars <- nCells
    } else {
        labelPars2 <- lapply(labelPars, "gparToList")
        nLabelPars <- length(labelPars)
    }
    parToLabelRatio <- nLabelPars / nLabels
    if (nLabelPars > nLabels && round(parToLabelRatio) == parToLabelRatio) {
        labels2 <- rep(labels2, length.out = nLabelPars)
        labelsNA <- rep_len(labelsNA, nLabelPars)
        nLabels <- nLabelPars
    } else if (nLabelPars != nLabels) {
        labels2 <- rep(labels2, length.out = nCells)
        labelsNA <- rep_len(labelsNA, nCells)
        labelPars2 <- rep_len(labelPars2, nCells)
        nLabels <- nCells
    }
    fsPar <- unname(unlist(get.gpar("fontsize")))
    maxLabelFS <- fsPar
    for (i in seq_along(labelPars2)) {
        thisSize <- labelPars2[[i]][["fontsize"]]
        if (is.null(thisSize)) {
            labelPars2[[i]][["fontsize"]] <- fsPar
        } else {
            maxLabelFS <- max(maxLabelFS, thisSize)
        }
    }
    axesParList <- gparToList(axesPar)
    if (!resizeAxes) {
        thisSize <- axesParList[["fontsize"]]
        if (is.null(thisSize)) {
            axesParList[["fontsize"]] <- fsPar
        }
        axesPar2 <- do.call("gpar", axesParList)
    }
    mainParList <- gparToList(mainPar)
    if (!resizeMain) {
        thisSize <- mainParList[["fontsize"]]
        if (is.null(thisSize)) {
            mainParList[["fontsize"]] <- fsPar
        }
        mainPar2 <- do.call("gpar", mainParList)
    }
    labParList <- gparToList(labPar)
    if (!resizeLab) {
        thisSize <- labParList[["fontsize"]]
        if (is.null(thisSize)) {
            labParList[["fontsize"]] <- fsPar
        }
        labPar2 <- do.call("gpar", labParList)
    }
    have.main <- length(main) > 0 && !all(vapply(main, allNA, TRUE,
                                                 USE.NAMES = FALSE))
    if (have.main) {
        if (is.expression(main)) {
            main2 <- list(main[1])
        } else {
            main2 <- as.list(main)[1]
        }
    }
    have.xlab <- length(xlab) > 0 && !all(vapply(xlab, allNA, TRUE,
                                                 USE.NAMES = FALSE))
    if (have.xlab) {
        if (is.expression(xlab)) {
            xlab2 <- list(xlab[1])
        } else {
            xlab2 <- as.list(xlab)[1]
        }
    }
    have.ylab <- length(ylab) > 0 && !all(vapply(ylab, allNA, TRUE,
                                                 USE.NAMES = FALSE))
    if (have.ylab) {
        if (is.expression(ylab)) {
            ylab2 <- list(ylab[1])
        } else {
            ylab2 <- as.list(ylab)[1]
        }
    }

    ## Argument checking
    if (have.xAxisLabels || have.yAxisLabels) {
        if (inherits(axesMargin, "unit")) {
            if (resizeAxes) {
                stop(gettextf("when '%s' is TRUE, '%s' must be a non-negative number",
                              "resizeAxes", "axesMargin", domain = "R-sisal"),
                     domain = NA)
            }
        } else if (!(is.numeric(axesMargin) && length(axesMargin) == 1 &&
                     is.finite(axesMargin) && axesMargin >= 0)) {
            stop(gettextf("if '%s' is not a unit, it must be a non-negative number",
                          "axesMargin", domain = "R-sisal"),
                 domain = NA)
        }
    }
    if (have.xlab) {
        if (inherits(xlabMargin, "unit")) {
            if (resizeLab) {
                stop(gettextf("when '%s' is TRUE, '%s' must be a non-negative number",
                              "resizeLab", "xlabMargin", domain = "R-sisal"),
                     domain = NA)
            }
        } else if (!is.numeric(xlabMargin) || length(xlabMargin) != 1 ||
                   !is.finite(xlabMargin) || xlabMargin < 0) {
            stop(gettextf("if '%s' is not a unit, it must be a non-negative number",
                          "xlabMargin", domain = "R-sisal"),
                 domain = NA)
        }
    }
    if (have.ylab) {
        if (inherits(ylabMargin, "unit")) {
            if (resizeLab) {
                stop(gettextf("when '%s' is TRUE, '%s' must be a non-negative number",
                              "resizeLab", "ylabMargin", domain = "R-sisal"),
                     domain = NA)
            }
        } else if (!is.numeric(ylabMargin) || length(ylabMargin) != 1 ||
                   !is.finite(ylabMargin) || ylabMargin < 0) {
            stop(gettextf("if '%s' is not a unit, it must be a non-negative number",
                          "ylabMargin", domain = "R-sisal"),
                 domain = NA)
        }
    }
    if (have.main) {
        if (inherits(mainMargin, "unit")) {
            if (resizeMain) {
                stop(gettextf("when '%s' is TRUE, '%s' must be a non-negative number",
                              "resizeMain", "mainMargin", domain = "R-sisal"),
                     domain = NA)
            }
        } else if (!is.numeric(mainMargin) || length(mainMargin) != 1 ||
                   !is.finite(mainMargin) || mainMargin < 0) {
            stop(gettextf("if '%s' is not a unit, it must be a non-negative number",
                          "mainMargin", domain = "R-sisal"),
                 domain = NA)
        }
    }

    ## Compute the actual sizes of different text elements and the
    ## relative sizes of different elements in the plot, taking into
    ## account some font size parameters etc.
    xy <- tableMeasures(xAxisLabels2[!xAxisNA], yAxisLabels2[!yAxisNA],
                        labels2[!labelsNA], targetAxesRatio = axesSize,
                        gpAxes = axesParList, gpLabels = labelPars2[!labelsNA],
                        main = if (have.main) main2,
                        xlab = if (have.xlab) xlab2,
                        ylab = if (have.ylab) ylab2,
                        gpMain = mainParList, gpLab = labParList,
                        nRows = nRows, nCols = nCols,
                        xlabExpand = xlabSize, ylabExpand = ylabSize,
                        mainExpand = mainSize, scaleLabels = forceAxesSize,
                        xAxisRot = xAxisRot, yAxisRot = yAxisRot,
                        xAxisJust = xAxisJust, yAxisJust = yAxisJust)
    labelJust <- rep.int(0.5, nLabels)
    labelsNA[!labelsNA] <- !xy[["labelNZ"]]
    labelsNotNA <- !labelsNA
    have.labels <- any(labelsNotNA)
    xAxisNA[!xAxisNA] <- !xy[["xAxisNZ"]]
    have.xAxisLabels <- !all(xAxisNA)
    yAxisNA[!yAxisNA] <- !xy[["yAxisNZ"]]
    have.yAxisLabels <- !all(yAxisNA)
    have.main <- xy[["mainNZ"]]
    have.xlab <- xy[["xlabNZ"]]
    have.ylab <- xy[["ylabNZ"]]
    internalLabels <- c(labels2[labelsNotNA], "x")
    internalLabelPars <- c(labelPars2[labelsNotNA],
                           list(list(fontsize = maxLabelFS)))
    if (have.labels && (!resizeLabels || !resizeTable)) {
        justTmp <- dynTextGrob(internalLabels, adjustJust = TRUE,
                               gp = internalLabelPars, resize = FALSE,
                               takeMeasurements = TRUE)[["vjust"]]
        labelJust[labelsNotNA] <- justTmp[1:(length(justTmp) - 1)]
    }

    ## Size of xlab, ylab margin (fixed size)
    resizeXlabMargin <- TRUE
    resizeYlabMargin <- TRUE
    if (!resizeLab) {
        if (have.xlab) {
            mfs <- dynTextGrob(xlab2, gp = labPar2,
                               resize = FALSE, takeMeasurements = TRUE)
            xlabWidthInches <- mfs[["sizingWidth"]]
            xlabHeightInches <- mfs[["sizingHeight"]]
            xlabHeight <- unit(xlabHeightInches, "inches")
        }
        if (have.ylab) {
            mfs <- dynTextGrob(ylab2, rot = 90, gp = labPar2,
                               resize = FALSE, takeMeasurements = TRUE)
            ylabWidthInches <- mfs[["sizingWidth"]]
            ylabHeightInches <- mfs[["sizingHeight"]]
            ylabWidth <- unit(ylabWidthInches, "inches")
        }
        if (inherits(xlabMargin, "unit")) {
            resizeXlabMargin <- FALSE
            xlabMarInches <-
                convertUnit(xlabMargin, unitTo = "inches",
                            axisFrom = "x", typeFrom = "dimension",
                            valueOnly = TRUE)
            xlabMar <- unit(xlabMarInches, "inches")
        }
        if (inherits(ylabMargin, "unit")) {
            resizeYlabMargin <- FALSE
            ylabMarInches <-
                convertUnit(ylabMargin, unitTo = "inches",
                            axisFrom = "y", typeFrom = "dimension",
                            valueOnly = TRUE)
            ylabMar <- unit(ylabMarInches, "inches")
        }
    }

    ## Size of x and y axes and their margin (fixed size)
    resizeAxesMargin <- TRUE
    xAxisNotNA <- !xAxisNA
    yAxisNotNA <- !yAxisNA
    if (!resizeAxes) {
        maxXAxisHeightInches <- 0
        xAxisLabelLeft <- numeric(nCols)
        xAxisLabelRight <- numeric(nCols)
        xAxisHJust <- numeric(nCols)
        xAxisVJust <- numeric(nCols)
        xAxisXShift <- numeric(nCols)
        xAxisYShift <- numeric(nCols)
        if (have.xAxisLabels) {
            xMeasurements <-
                dynTextGrob(xAxisLabels[xAxisNotNA], adjustJust = TRUE,
                            gp = axesPar2, resize = FALSE, rotJust = TRUE,
                            rot = xAxisRot, just = xAxisJust,
                            takeMeasurements = TRUE)
            maxXAxisHeightInches <- xMeasurements[["sizingHeight"]]
            xAxisLabelLeft[xAxisNotNA] <- xMeasurements[["left"]]
            xAxisLabelRight[xAxisNotNA] <- xMeasurements[["right"]]
            xAxisHJust[xAxisNotNA] <- xMeasurements[["hjust"]]
            xAxisVJust[xAxisNotNA] <- xMeasurements[["vjust"]]
            xAxisXShift[xAxisNotNA] <- xMeasurements[["xShift"]]
            xAxisYShift[xAxisNotNA] <- xMeasurements[["yShift"]]
            xAxisX2 <- xAxisX + xAxisXShift
            xAxisY2 <- xAxisY + xAxisYShift
        }
        maxYAxisWidthInches <- 0
        yAxisLabelBottom <- numeric(nRows)
        yAxisLabelTop <- numeric(nRows)
        yAxisHJust <- numeric(nRows)
        yAxisVJust <- numeric(nRows)
        yAxisXShift <- numeric(nRows)
        yAxisYShift <- numeric(nRows)
        if (have.yAxisLabels) {
            yMeasurements <-
                dynTextGrob(yAxisLabels[yAxisNotNA], adjustJust = TRUE,
                            gp = axesPar2, resize = FALSE, rotJust = TRUE,
                            rot = yAxisRot, just = yAxisJust,
                            takeMeasurements = TRUE)
            maxYAxisWidthInches <- yMeasurements[["sizingWidth"]]
            yAxisLabelBottom[yAxisNotNA] <- yMeasurements[["bottom"]]
            yAxisLabelTop[yAxisNotNA] <- yMeasurements[["top"]]
            yAxisHJust[yAxisNotNA] <- yMeasurements[["hjust"]]
            yAxisVJust[yAxisNotNA] <- yMeasurements[["vjust"]]
            yAxisXShift[yAxisNotNA] <- yMeasurements[["xShift"]]
            yAxisYShift[yAxisNotNA] <- yMeasurements[["yShift"]]
            yAxisX2 <- yAxisX + yAxisXShift
            yAxisY2 <- yAxisY + yAxisYShift
        }
        if (have.xAxisLabels) {
            xAxisSide <- unit(maxXAxisHeightInches, "inches")
        }
        if (have.yAxisLabels) {
            yAxisSide <- unit(maxYAxisWidthInches, "inches")
        }
        if (inherits(axesMargin, "unit")) {
            resizeAxesMargin <- FALSE
            xyMarginSideInches <- minConvertXY(axesMargin)
            xyMarginSide <- unit(xyMarginSideInches, "inches")
        }
    }

    ## Size of main title and its margin (fixed size)
    resizeMainMargin <- TRUE
    if (have.main && !resizeMain) {
        mfs <- dynTextGrob(main2, gp = mainPar2,
                           resize = FALSE, takeMeasurements = TRUE)
        mainWidthInches <- mfs[["sizingWidth"]]
        mainHeightInches <- mfs[["sizingHeight"]]
        mainHeight <- unit(mainHeightInches, "inches")
        if (inherits(mainMargin, "unit")) {
            resizeMainMargin <- FALSE
            mainMarginInches <- minConvertXY(mainMargin)
            mainMargin2 <- unit(mainMarginInches, "inches")
        }
    }

    ## Layout, number of rows and columns
    layRows <- 1 + 2 * have.main + nRows + 2 * have.xAxisLabels +
        2 * have.xlab
    layCols <- 1 + 2 * have.ylab + 2 * have.yAxisLabels + nCols

    ## Sums of width and height of fixed-size ("rigid") objects
    rigidWidthInches <- 0
    rigidHeightInches <- 0
    ## Width and height of the table when resizeTable is FALSE (adds
    ## to rigidWidthInches and rigidHeightInches)
    if (!resizeTable) {
        labelsCex <- FALSE
        if (have.labels) {
            whichLabels <- which(labelsNotNA)
            nNotNA <- length(whichLabels)
            gpS <- labelPars2[whichLabels]
            theCex <- numeric(nNotNA)
            for (i in seq_len(nNotNA)) {
                thisCex <- gpS[[i]][["cex"]]
                if (is.numeric(thisCex) && length(thisCex) == 1) {
                    theCex[i] <- thisCex
                }
                gpS[[i]][["cex"]] <- NULL
            }
            mfs <- dynTextGrob(labels2[whichLabels], takeMeasurements = TRUE,
                               gp = gpS, resize = FALSE, adjustJust = FALSE,
                               vjust = labelJust[whichLabels])
            cellSideInches <- max(unlist(mfs[c("sizingWidth", "sizingHeight")]))
            if (!anyResize) {
                labelHalfWidths <- numeric(length(labelsNA))
                labelTops <- numeric(length(labelsNA))
                labelBottoms <- numeric(length(labelsNA))
                bigCex <- logical(length(labelsNA))
                labelHalfWidths[whichLabels] <- mfs[["right"]]
                labelTops[whichLabels] <- mfs[["top"]]
                labelBottoms[whichLabels] <- mfs[["bottom"]]
                idx <- which(theCex > 1)
                labelsCex <- length(idx) > 0
                if (labelsCex) {
                    bigCex[idx] <- TRUE
                    theCex <- theCex[idx]
                    labelHalfWidths[idx] <- theCex * labelHalfWidths[idx]
                    labelTops[idx] <- theCex * labelTops[idx]
                    labelBottoms[idx] <- theCex * labelBottoms[idx]
                }
            }
        } else {
            tmpGrob <- textGrob("x", gp = gpar(fontsize = maxLabelFS))
            cellSideInches <-
                max(convertWidth(grobWidth(tmpGrob),
                                 unitTo = "inches", valueOnly = TRUE),
                    convertHeight(grobHeight(tmpGrob),
                                  unitTo = "inches", valueOnly = TRUE))
        }

        cellSide <- unit(cellSideInches, "inches")
        xCell <- nCols * cellSideInches
        yCell <- nRows * cellSideInches
        rigidWidthInches <- rigidWidthInches + xCell
        rigidHeightInches <- rigidHeightInches + yCell
        if (!resizeAxes && resizeAxesMargin) {
            resizeAxesMargin <- FALSE
            xyMarginSideInches <- axesMargin * cellSideInches
            xyMarginSide <- unit(xyMarginSideInches, "inches")
        }
        if (have.main && !resizeMain && resizeMainMargin) {
            resizeMainMargin <- FALSE
            mainMarginInches <- mainMargin * cellSideInches
            mainMargin2 <- unit(mainMarginInches, "inches")
        }
        if (!resizeLab && resizeXlabMargin) {
            resizeXlabMargin <- FALSE
            xlabMarInches <- xlabMargin * cellSideInches
            xlabMar <- unit(xlabMarInches, "inches")
        }
        if (!resizeLab && resizeYlabMargin) {
            resizeYlabMargin <- FALSE
            ylabMarInches <- ylabMargin * cellSideInches
            ylabMar <- unit(ylabMarInches, "inches")
        }
    }
    if (have.ylab && !resizeLab) {
        rigidWidthInches <- rigidWidthInches + ylabWidthInches
        if (!resizeYlabMargin) {
            rigidWidthInches <- rigidWidthInches + ylabMarInches
        }
    }
    if (have.yAxisLabels && !resizeAxes) {
        rigidWidthInches <- rigidWidthInches + maxYAxisWidthInches
        if (!resizeAxesMargin) {
            rigidWidthInches <- rigidWidthInches + xyMarginSideInches
        }
    }
    if (have.main && !resizeMain) {
        rigidHeightInches <- rigidHeightInches + mainHeightInches
        if (!resizeMainMargin) {
            rigidHeightInches <- rigidHeightInches + mainMarginInches
        }
    }
    if (have.xAxisLabels && !resizeAxes) {
        rigidHeightInches <- rigidHeightInches + maxXAxisHeightInches
        if (!resizeAxesMargin) {
            rigidHeightInches <- rigidHeightInches + xyMarginSideInches
        }
    }
    if (have.xlab && !resizeLab) {
        rigidHeightInches <- rigidHeightInches + xlabHeightInches
        if (!resizeXlabMargin) {
            rigidHeightInches <- rigidHeightInches + xlabMarInches
        }
    }

    ## Sizes of dummGrob and objects depending on it can be constant or dynamic
    if (inherits(width, "unit") && inherits(height, "unit")) {
        hSize <- width
        vSize <- height
    } else {
        hSize <- unit(width, default.units)
        vSize <- unit(height, default.units)
    }
    ## This is never drawn
    dummGrob <- rectGrob(width = max(unit(0, "inches"),
                             hSize - unit(rigidWidthInches, "inches")),
                         height = max(unit(0, "inches"),
                             vSize - unit(rigidHeightInches, "inches")))

    ## hUnits and vUnits are the total horizontal and vertical size of
    ## table cells whose size depends on the size of dummGrob, where 1
    ## unit is the size of a square cell in the nRows * nCols contents
    ## section (if resizeTable is TRUE)
    if (resizeTable) {
        hUnits <- nCols
    } else {
        hUnits <- 0
    }
    if (have.yAxisLabels) {
        if (resizeAxes) {
            hUnits <- hUnits + xy[["yAxisRatio"]]
        }
        if (resizeAxesMargin) {
            hUnits <- hUnits + axesMargin
        }
    }
    if (have.ylab) {
        if (resizeLab) {
            hUnits <- hUnits + xy[["ylabRatio"]]
        }
        if (resizeYlabMargin) {
            hUnits <- hUnits + ylabMargin
        }
    }
    if (resizeTable) {
        vUnits <- nRows
    } else {
        vUnits <- 0
    }
    if (have.xAxisLabels) {
        if (resizeAxes) {
            vUnits <- vUnits + xy[["xAxisRatio"]]
        }
        if (resizeAxesMargin) {
            vUnits <- vUnits + axesMargin
        }
    }
    if (have.xlab) {
        if (resizeLab) {
            vUnits <- vUnits + xy[["xlabRatio"]]
        }
        if (resizeXlabMargin) {
            vUnits <- vUnits + xlabMargin
        }
    }
    if (have.main) {
        if (resizeMain) {
            vUnits <- vUnits + xy[["mainRatio"]]
        }
        if (resizeMainMargin) {
            vUnits <- vUnits + mainMargin
        }
    }

    ## Width and / or height of graphical objects as units
    if (resizeTable) {
        cellSide <- min(unit(1 / hUnits, "grobwidth", data = dummGrob),
                        unit(1 / vUnits, "grobheight", data = dummGrob))
    }
    xAxisWidth <- xy[["xAxisWidth"]]
    xAxisHeight <- xy[["xAxisHeight"]]
    yAxisWidth <- xy[["yAxisWidth"]]
    yAxisHeight <- xy[["yAxisHeight"]]
    if (resizeAxes) {
        xyMarginSide <- axesMargin * cellSide
        if (have.xAxisLabels) {
            xAxisSide <-
                min(xAxisHeight / xAxisWidth * cellSide,
                    unit(xy[["xAxisRatio"]] / hUnits, "grobwidth",
                         data = dummGrob),
                    unit(xy[["xAxisRatio"]] / vUnits, "grobheight",
                         data = dummGrob))
        }
        if (have.yAxisLabels) {
            yAxisSide <-
                min(yAxisWidth / yAxisHeight * cellSide,
                    unit(xy[["yAxisRatio"]] / hUnits, "grobwidth",
                         data = dummGrob),
                    unit(xy[["yAxisRatio"]] / vUnits, "grobheight",
                         data = dummGrob))
        }
    }
    layHeights <- vector(mode = "list", length = layRows)
    layHeights[[1]] <- unit(0, "inches")
    rowCount <- 1
    if (have.main) {
        if (resizeMain) {
            mainMargin2 <- mainMargin * cellSide
            mainHeight <- min(xy[["mainHeight"]] / xy[["mainWidth"]] *
                              nCols * cellSide,
                              unit(xy[["mainRatio"]] / hUnits, "grobwidth",
                                   data = dummGrob),
                              unit(xy[["mainRatio"]] / vUnits, "grobheight",
                                   data = dummGrob))
        }
        layHeights[rowCount + 1:2] <- list(mainHeight, mainMargin2)
        rowCount <- rowCount + 2
    }
    layHeights[rowCount + seq_len(nRows)] <- rep(list(cellSide), nRows)
    rowCount <- rowCount + nRows
    if (have.xAxisLabels) {
        layHeights[rowCount + 1:2] <- list(xyMarginSide, xAxisSide)
        rowCount <- rowCount + 2
    } else {
        xAxisSide <- NULL
    }
    if (resizeLab) {
        if (have.xlab) {
            xlabMar <- xlabMargin * cellSide
            xlabHeight <- min(xy[["xlabHeight"]] / xy[["xlabWidth"]] *
                              nCols * cellSide,
                              unit(xy[["xlabRatio"]] / hUnits, "grobwidth",
                                   data = dummGrob),
                              unit(xy[["xlabRatio"]] / vUnits, "grobheight",
                                   data = dummGrob))
        }
        if (have.ylab) {
            ylabMar <- ylabMargin * cellSide
            ylabWidth <- min(xy[["ylabWidth"]] / xy[["ylabHeight"]] *
                             nRows * cellSide,
                             unit(xy[["ylabRatio"]] / hUnits, "grobwidth",
                                  data = dummGrob),
                             unit(xy[["ylabRatio"]] / vUnits, "grobheight",
                                  data = dummGrob))
        }
    }
    if (have.xlab) {
        layHeights[rowCount + 1:2] <- list(xlabMar, xlabHeight)
    }
    layHeights <- do.call("unit.c", layHeights)
    layWidths <- vector(mode = "list", length = layCols)
    layWidths[[1]] <- unit(0, "inches")
    colCount <- 1
    if (have.ylab) {
        layWidths[colCount + 1:2] <- list(ylabWidth, ylabMar)
        colCount <- colCount + 2
    }
    if (have.yAxisLabels) {
        layWidths[colCount + 1:2] <- list(yAxisSide, xyMarginSide)
        colCount <- colCount + 2
    } else {
        yAxisSide <- NULL
    }
    layWidths[colCount + seq_len(nCols)] <- rep(list(cellSide), nCols)
    layWidths <- do.call("unit.c", layWidths)

    ## Lists for grobs and viewports, with some excess length.
    ## The drawing order will be innerGrobs before outerGrobs.
    outerGrobs <- vector(mode = "list",
                         length = layRows * layCols - nRows * nCols)
    innerGrobs <- vector(mode = "list",
                         length = nRows * nCols + outerRect +
                         innerLines * ((nRows-1) * nCols + (nCols-1) * nRows))
    vps <- vector(mode = "list", length = layRows * layCols + outerRect)
    ## Counts how many grobs and viewports have been added to the lists
    vpCounter <- 0
    outerCount <- 0
    innerCount <- 0

### Create grobs and viewports, add to the lists
    if (any(bgNA) || any(!stripeNA[!bgNA])) {
        stripeWidth2 <- stripeWidth / nStripes
        stripeSep <- 1 / nStripes
        stripeY <- seq(from=stripeSep / 2, by = stripeSep, length.out=nStripes)
        stripeGrobNA <-
            grobTree(rectGrob(gp = gpar(fill=naFill, col=naFill),
                              name = "stripeBg"),
                     stripePolygon(y = stripeY, rot = stripeRot,
                                   width = stripeWidth2,
                                   xyScale = stripeScale,
                                   default.units = "npc",
                                   name = "stripes"),
                     gp = gpar(col = naStripes, fill = naStripes))
    }
    xAxisGrobs <- vector(mode = "list", length = nXAxis)
    yAxisGrobs <- vector(mode = "list", length = nYAxis)
    if (resizeAxes) {
        for (i in which(xAxisNotNA)) {
            xAxisGrobs[[i]] <-
                textGrob(xAxisLabels2[[i]], gp = axesPar,
                         rot = xAxisRot, x = xAxisX, y = xAxisY)
        }
        for (i in which(yAxisNotNA)) {
            yAxisGrobs[[i]] <-
                textGrob(yAxisLabels2[[i]], gp = axesPar,
                         rot = yAxisRot, x = yAxisX, y = yAxisY)
        }
    } else {
        for (i in which(xAxisNotNA)) {
            xAxisGrobs[[i]] <-
                textGrob(xAxisLabels2[[i]], gp = axesPar2,
                         rot = xAxisRot, x = xAxisX2[i], y = xAxisY2[i],
                         hjust = xAxisHJust[i], vjust = xAxisVJust[i])
        }
        for (i in which(yAxisNotNA)) {
            yAxisGrobs[[i]] <-
                textGrob(yAxisLabels2[[i]], gp = axesPar2,
                         rot = yAxisRot, x = yAxisX2[i], y = yAxisY2[i],
                         hjust = yAxisHJust[i], vjust = yAxisVJust[i])
        }
    }

    ## Row by row...
    cellNames <- twoPartName(rep.int(seq_len(nRows), nCols),
                             rep(seq_len(nCols), each = nRows))
    bgNames <- twoPartNameChar("bg", cellNames)
    textNames <- twoPartNameChar("text", cellNames)
    treeNames <- twoPartNameChar("tree", cellNames)
    yNames <- twoPartName(seq_len(nRows), "yaxis")
    for (Row in seq_len(nRows)) {
        Row2 <- 1 + Row + 2 * have.main
        ## Y axis label
        yAxisMod <- (Row - 1) %% nYAxis + 1
        if (yAxisNotNA[yAxisMod]) {
            vpCounter <- vpCounter + 1
            outerCount <- outerCount + 1
            Col2 <- 2 + 2 * have.ylab
            thisName <- yNames[Row]
            vps[[vpCounter]] <-
                viewport(layout.pos.row = Row2, layout.pos.col = Col2,
                         name = thisName)
            outerGrobs[[outerCount]] <-
                editGrob(yAxisGrobs[[yAxisMod]], name = thisName,
                         vp = vpPath("Parent", thisName))
        }
        Col2 <- 2 + 2 * have.yAxisLabels + 2 * have.ylab
        ## Actual contents of the table
        for (Col in seq_len(nCols)) {
            linearRowCol <- Row - 1 + (Col - 1) * nRows
            vpCounter <- vpCounter + 1
            innerCount <- innerCount + 1
            thisName <- cellNames[linearRowCol + 1]
            vps[[vpCounter]] <-
                viewport(layout.pos.row = Row2, layout.pos.col = Col2,
                         name = thisName)
            vpP <- vpPath("Parent", thisName)
            labelsMod <- linearRowCol %% nLabels + 1
            labelPar <- labelPars2[[labelsMod]]
            labelPar[["col"]] <- fg2[[linearRowCol %% nFg + 1]]
            bgMod <- linearRowCol + 1
            if (bgNA[bgMod]) {
                bgGrob <- editGrob(stripeGrobNA, vp = vpP,
                                   name = bgNames[linearRowCol + 1])
            } else if (stripeNA[bgMod]) {
                thisBg <- bg2[[bgMod]]
                bgGrob <- rectGrob(gp = gpar(fill = thisBg, col = thisBg),
                                   vp = vpP,
                                   name = bgNames[linearRowCol + 1])
            } else {
                thisBg1 <- bg2[[bgMod]]
                thisBg2 <- stripeCol2[[bgMod]]
                bgGrob <- editGrob(stripeGrobNA, vp = vpP,
                                   name = bgNames[linearRowCol + 1],
                                   gp = gpar(col = thisBg2, fill = thisBg2))
                bgGrob <- editGrob(bgGrob, gPath = "stripeBg",
                                   gp = gpar(col = thisBg1, fill = thisBg1))
            }
            if (resizeLabels) {
                labelPar["fontsize"] <- NULL
            }
            if (labelsNA[labelsMod]) {
                fgGrob <- NULL
            } else {
                fgGrob <- textGrob(labels2[[labelsMod]],
                                   gp = do.call("gpar", labelPar), vp = vpP,
                                   name = textNames[linearRowCol + 1],
                                   vjust = labelJust[labelsMod])
            }
            innerGrobs[[innerCount]] <-
                grobTree(bgGrob, fgGrob, name = treeNames[linearRowCol + 1])
            Col2 <- Col2 + 1
        }
    }
    ## X axis labels
    Row2 <- nRows + 3 + 2 * have.main
    Col2 <- 1 + 2 * have.yAxisLabels + 2 * have.ylab
    xNames <- twoPartName("xaxis", seq_len(nCols))
    for (Col in seq_len(nCols)[xAxisNotNA]) {
        vpCounter <- vpCounter + 1
        outerCount <- outerCount + 1
        thisName <- xNames[Col]
        vps[[vpCounter]] <- viewport(layout.pos.row = Row2, name = thisName,
                                     layout.pos.col = Col2 + Col)
        outerGrobs[[outerCount]] <-
            editGrob(xAxisGrobs[[(Col - 1) %% nXAxis + 1]], name = thisName,
                     vp = vpPath("Parent", thisName))
    }
    ## X label
    if (have.xlab) {
        vpCounter <- vpCounter + 1
        outerCount <- outerCount + 1
        xlabRow <- nRows + 2 * have.xAxisLabels + 3 + 2 * have.main
        firstCol <- 2 + 2 * have.yAxisLabels + 2 * have.ylab
        lastCol <- firstCol + nCols - 1
        xlabCol <- c(firstCol, lastCol)
        vps[[vpCounter]] <-
            viewport(layout.pos.row = xlabRow, layout.pos.col = xlabCol,
                     name = "xlab")
        if (resizeLab) {
            outerGrobs[[outerCount]] <-
                dynTextGrob(xlab2[[1]], gp = labPar,
                            sizingWidth = xy[["xlabWidth"]],
                            sizingHeight = xy[["xlabHeight"]],
                            name = "xlab", vp = vpPath("Parent", "xlab"),)
        } else {
            outerGrobs[[outerCount]] <-
                textGrob(xlab2[[1]], gp = labPar2, name = "xlab",
                         vp = vpPath("Parent", "xlab"))
        }
    }
    ## Y label
    if (have.ylab) {
        vpCounter <- vpCounter + 1
        outerCount <- outerCount + 1
        ylabCol <- 2
        firstRow <- 2 + 2 * have.main
        lastRow <- firstRow + nRows - 1
        ylabRow <- c(firstRow, lastRow)
        vps[[vpCounter]] <-
            viewport(layout.pos.row = ylabRow, layout.pos.col = ylabCol,
                     name = "ylab")
        if (resizeLab) {
            outerGrobs[[outerCount]] <-
                dynTextGrob(ylab2[[1]], gp = labPar, rot = 90, name = "ylab",
                            vp = vpPath("Parent", "ylab"),
                            sizingWidth = xy[["ylabWidth"]],
                            sizingHeight = xy[["ylabHeight"]])
        } else {
            outerGrobs[[outerCount]] <-
                textGrob(ylab2[[1]], gp = labPar2, rot = 90, name = "ylab",
                         vp = vpPath("Parent", "ylab"))
        }
    }
    ## Main title
    if (have.main) {
        vpCounter <- vpCounter + 1
        outerCount <- outerCount + 1
        firstCol <- 2 + 2 * have.yAxisLabels + 2 * have.ylab
        lastCol <- firstCol + nCols - 1
        mainCol <- c(firstCol, lastCol)
        mainRow <- 2
        vps[[vpCounter]] <-
            viewport(layout.pos.row = mainRow, layout.pos.col = mainCol,
                     name = "main")
        if (resizeMain) {
            outerGrobs[[outerCount]] <-
                dynTextGrob(main2[[1]],
                            gp = mainPar, vp = vpPath("Parent", "main"),
                            sizingWidth = xy[["mainWidth"]],
                            sizingHeight = xy[["mainHeight"]],
                            name = "main")
        } else {
            outerGrobs[[outerCount]] <- textGrob(main2[[1]], gp = mainPar2,
                                                 vp = vpPath("Parent", "main"),
                                                 name = "main")
        }
    }

    ## Line segments between table cells: Maximize contrast to
    ## neighboring cells (black, white and different shades of gray
    ## available)
    if (innerLines) {
        linePar2 <- gparToList(linePar)
        colorMissing <- !("col" %in% names(linePar2))
        if (colorMissing) {
            linePar2[["col"]] <- "black"
            colorIdx <- which(names(linePar2) == "col")[1]
        }
        ## Horizontal
        hLineNames <- twoPartNameChar("hLine", cellNames)
        for (Row in seq_len(nRows - 1)) {
            for (Col in seq_len(nCols)) {
                linearRowCol <- c(Row, Row + 1) - 1 + (Col - 1) * nRows + 1
                innerCount <- innerCount + 1
                vpName <- cellNames[linearRowCol[1]]
                if (colorMissing) {
                    linePar2[[colorIdx]] <- fgGrey(bgLuma[linearRowCol])
                }
                innerGrobs[[innerCount]] <-
                    segmentsGrob(x0 = 0, y0 = 0, x1 = 1, y1 = 0,
                                 gp = do.call("gpar", linePar2),
                                 name = hLineNames[linearRowCol[1]],
                                 vp = vpPath("Parent", vpName))
            }
        }
        ## Vertical
        vLineNames <- twoPartNameChar("vLine", cellNames)
        for (Col in seq_len(nCols - 1)) {
            for (Row in seq_len(nRows)) {
                linearRowCol <- Row - 1 + (c(Col, Col + 1) - 1) * nRows + 1
                innerCount <- innerCount + 1
                vpName <- cellNames[linearRowCol[1]]
                if (colorMissing) {
                    linePar2[[colorIdx]] <- fgGrey(bgLuma[linearRowCol])
                }
                innerGrobs[[innerCount]] <-
                    segmentsGrob(x0 = 1, y0 = 0, x1 = 1, y1 = 1,
                                 gp = do.call("gpar", linePar2),
                                 name = vLineNames[linearRowCol[1]],
                                 vp = vpPath("Parent", vpName))
            }
        }
    }

    ## Frame around table
    if (outerRect) {
        rectPar2 <- gparToList(rectPar)
        rectPar2[["fill"]] <- 0 # transparent
        vpCounter <- vpCounter + 1
        innerCount <- innerCount + 1
        Row <- seq(from = 2 + 2 * have.main, by = 1, length.out = nRows)
        Col <- seq(from = 2 + 2 * have.yAxisLabels + 2 * have.ylab,
                   by = 1, length.out = nCols)
        vps[[vpCounter]] <- viewport(layout.pos.row = Row,
                                     layout.pos.col = Col, name = "outerRect")
        innerGrobs[[innerCount]] <- rectGrob(gp = do.call("gpar", rectPar2),
                                             name = "outerRect",
                                             vp = vpPath("Parent", "outerRect"))
    }

    ## Compute numeric 'hjust' and 'vjust'
    tmpJust <- justHV(just)
    hjust <- tmpJust[[1]]
    vjust <- tmpJust[[2]]
    ## Convert 'x' and 'y' to unit objects
    if (inherits(x, "unit")) {
        X <- x
    } else {
        X <- unit(x, default.units)
    }
    if (inherits(y, "unit")) {
        Y <- y
    } else {
        Y <- unit(y, default.units)
    }

    ## Create parent viewport
    if (anyResize) {
        vpWidth <- sum(layWidths)
        vpHeight <- sum(layHeights)
    } else {
        ## In this branch, everything has fixed size. Therefore it is
        ## possible to check if the text elements will overflow the
        ## layout and by how much. The margin is then widened to cover
        ## all the oversized elements..
        leftHang <- 0
        rightHang <- 0
        if ((have.main && mainWidthInches > xCell) ||
            (have.xlab && xlabWidthInches > xCell) ||
            have.xAxisLabels) {
            xLeft <- 0
            if (have.ylab) {
                xLeft <- xLeft + ylabWidthInches + ylabMarInches
            }
            if (have.yAxisLabels) {
                xLeft <- xLeft + maxYAxisWidthInches + xyMarginSideInches
            }
            xCenter <- xLeft + xCell / 2
            if (have.main) {
                leftHang <- max(leftHang, mainWidthInches / 2 - xCenter)
                rightHang <-
                    max(rightHang,
                        xCenter + mainWidthInches / 2 - rigidWidthInches)
            }
            if (have.xlab) {
                leftHang <- max(leftHang, xlabWidthInches / 2 - xCenter)
                rightHang <-
                    max(rightHang,
                        xCenter + xlabWidthInches / 2 - rigidWidthInches)
            }
            if (have.xAxisLabels) {
                colSeq <- cellSideInches *
                    (xAxisX2 + seq(from = 0, by = 1, length.out = nCols))
                xAxisLabelLeftMaxHang <-
                    xLeft - min(colSeq + xAxisLabelLeft)
                xAxisLabelRightMaxHang <-
                    max(colSeq + xAxisLabelRight) - xCell
                leftHang <- max(leftHang, xAxisLabelLeftMaxHang)
                rightHang <- max(rightHang, xAxisLabelRightMaxHang)
            }
        }
        bottomHang <- 0
        topHang <- 0
        if ((have.ylab && ylabHeightInches > yCell) || have.yAxisLabels) {
            yBottom <- 0
            if (have.xlab) {
                yBottom <- yBottom + xlabHeightInches + xlabMarInches
            }
            if (have.xAxisLabels) {
                yBottom <- yBottom + maxXAxisHeightInches + xyMarginSideInches
            }
            yCenter <- yBottom + yCell / 2
            if (have.ylab) {
                bottomHang <- max(bottomHang, ylabHeightInches / 2 - yCenter)
                topHang <-
                    max(topHang,
                        yCenter + ylabHeightInches / 2 - rigidHeightInches)
            }
            if (have.yAxisLabels) {
                rowSeq <- cellSideInches *
                    (yAxisY2 +
                     seq(from = nRows - 1, by = -1, length.out = nRows))
                yAxisLabelBottomMaxHang <-
                    yBottom - min(rowSeq + yAxisLabelBottom)
                yAxisLabelTopMaxHang <-
                    max(rowSeq + yAxisLabelTop) - yCell
                bottomHang <- max(bottomHang, yAxisLabelBottomMaxHang)
                topHang <- max(topHang, yAxisLabelTopMaxHang)
            }
            if (labelsCex) {
                bigHalfWidths <- labelHalfWidths[bigCex]
                bigTops <- labelTops[bigCex]
                bigBottoms <- labelBottoms[bigCex]
                bigCex2 <- rep(0, length(bigCex))
                nBigCex <- sum(bigCex)
                bigCex2[bigCex] <- 1:nBigCex
                bigCexLong <- rep_len(bigCex, nRows * nCols)
                bigCex2Long <- rep_len(bigCex2, nRows * nCols)[bigCexLong]
                rowNums <- rep(nRows:1, nCols)[bigCexLong]
                colNums <- rep(1:nCols, each = nRows)[bigCexLong]
                for (i in seq_len(nBigCex)) {
                    matches <- which(bigCex2Long == i)
                    rowMatches <- rowNums[matches]
                    colMatches <- colNums[matches]
                    halfWidth <- bigHalfWidths[i]
                    rowSeq <- (rowMatches - 0.5) * cellSideInches
                    colSeq <- (colMatches - 0.5) * cellSideInches
                    labelLeftMaxHang <- max(-colSeq) + halfWidth - xLeft
                    labelRightMaxHang <- max(colSeq) + halfWidth - xCell
                    labelBottomMaxHang <- yBottom - min(rowSeq) - bigBottoms[i]
                    labelTopMaxHang <- max(rowSeq) + bigTops[i] - yCell
                    leftHang <- max(leftHang, labelLeftMaxHang)
                    rightHang <- max(rightHang, labelRightMaxHang)
                    bottomHang <- max(bottomHang, labelBottomMaxHang)
                    topHang <- max(topHang, labelTopMaxHang)
                }
            }
        }
        vpWidth <- unit(rigidWidthInches + leftHang + rightHang, "inches")
        vpHeight <- unit(rigidHeightInches + bottomHang + topHang, "inches")
        if (leftHang > 0) {
            layWidths <- unit.c(unit(leftHang, "inches"), layWidths[-1])
        }
        if (rightHang > 0) {
            layCols <- layCols + 1
            layWidths <- unit.c(layWidths, unit(rightHang, "inches"))
        }
        if (topHang > 0) {
            layHeights <- unit.c(unit(topHang, "inches"), layHeights[-1])
        }
        if (bottomHang > 0) {
            layRows <- layRows + 1
            layHeights <- unit.c(layHeights, unit(bottomHang, "inches"))
        }
    }
    the.layout <- grid.layout(layRows, layCols,
                              widths = layWidths, heights = layHeights)
    parentVp <- viewport(x = X, y = Y, width = vpWidth, height = vpHeight,
                         just = just, clip = clip,
                         layout = the.layout, name = "Parent")
    ## A gList (the grobs) and a vpTree (the viewports) go to a gTree
    grobs <- do.call("gList", c(innerGrobs[seq_len(innerCount)],
                                outerGrobs[seq_len(outerCount)]))
    childVps <- do.call("vpList", vps[seq_len(vpCounter)])
    viewports <- vpTree(parent = parentVp, children = childVps)
    if (resizeLabels && resizeTable && have.labels) {
        nResizeLabels <- sum(labelsNotNA)
    } else {
        nResizeLabels <- 0
    }
    myPlot <-
        gTree(children = grobs, childrenvp = viewports,
              width = vpWidth, height = vpHeight, vpLayout = the.layout,
              x = X, y = Y, just = just, clip = clip,
              internals = list(hjust = hjust, vjust = vjust,
                  resizeAxes = resizeAxes, labelDim = xy[["labelDim"]],
                  xAxisWidth = xAxisWidth, xAxisHeight = xAxisHeight,
                  yAxisWidth = yAxisWidth, yAxisHeight = yAxisHeight,
                  cellSide = cellSide, xAxisSide = xAxisSide,
                  yAxisSide = yAxisSide, nRows = nRows, nCols = nCols,
                  labelsNotNA = labelsNotNA, xAxisNotNA = xAxisNotNA,
                  yAxisNotNA = yAxisNotNA, labels = internalLabels,
                  labelPars = internalLabelPars, axesPar = axesParList,
                  xAxisLabels = xAxisLabels2[xAxisNotNA],
                  yAxisLabels = yAxisLabels2[yAxisNotNA],
                  xAxisRot = xAxisRot, yAxisRot = yAxisRot,
                  xAxisJust = xAxisJust, yAxisJust = yAxisJust,
                  textNames = textNames, xNames = xNames, yNames = yNames,
                  nResizeLabels = nResizeLabels, childVps = childVps),
              cl = "sisalTable", name = name, gp = gp, vp = vp)
    if (draw) {
        grid.draw(myPlot)
    }
    ## An object of class c("sisalTable", "gTree", ...) is returned
    invisible(myPlot)
}
