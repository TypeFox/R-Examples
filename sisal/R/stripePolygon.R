### File R/stripePolygon.R
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

## Returns a grob describing stripes, by default reaching across the
## viewport.  Parameters x, y and width are value only, "npc" units by
## default (other units can be tried via default.units in ...).
## Rotation happens around (0.5, 0.5).  Resulting coordinates (in the
## 0..1 range) can be scaled, then shifted.
## All the arguments (except those in ...) can have more than one
## element. The arguments are recycled to the maximum length of all
## arguments (except xyScale and ...), which is taken to be the number
## of stripes.
stripePolygon <- function(x = NULL, y = NULL, width = 0.1, rot = -45,
                          xyScale = 1, xScale = xyScale, yScale = xyScale,
                          xShift = 0.5 - 0.5 * xScale,
                          yShift = 0.5 - 0.5 * yScale, ...) {
    nStripes <- max(length(width), length(rot), length(xScale),
                    length(yScale), length(xShift), length(yShift))
    if (!is.null(x)) {        # x takes precedence
        useX <- TRUE
        nStripes <- max(nStripes, length(x))
        x2 <- rep_len(x, nStripes)
    } else if (!is.null(y)) {
        useX <- FALSE
        nStripes <- max(nStripes, length(y))
        y2 <- rep_len(y, nStripes)
    } else {
        ## no coordinates given
        return(grob(...))     # "empty" grob
    }
    width2 <- rep_len(width, nStripes)
    xScale2 <- rep_len(xScale, nStripes)
    yScale2 <- rep_len(yScale, nStripes)
    xShift2 <- rep_len(xShift, nStripes)
    yShift2 <- rep_len(yShift, nStripes)
    rot2 <- round(rot) %% 360
    isPerpendicular <- rep_len(rot2 %% 90 == 0, nStripes)
    rot2 <- rot2 / 180 * pi
    sinRot <- rep_len(sin(rot2), nStripes)
    cosRot <- rep_len(cos(rot2), nStripes)
    ## x1, y1 are coordinates before rotation and clipping
    if (useX) {        # x takes precedence
        xLeft  <- pmin(0.5, pmax(-0.5, x2 - 0.5 * (width2 + 1)))
        xRight <- pmin(0.5, pmax(-0.5, x2 + 0.5 * (width2 - 1)))
        minX <- min(xLeft)
        maxX <- max(xRight)
        if (maxX == -0.5 || minX == 0.5) {
            ## outside drawing area
            return(grob(...)) # "empty" grob
        }
        y1 <- rep(c(-0.5, 0.5), each=2)
    } else {
        yBottom <- pmin(0.5, pmax(-0.5, y2 - 0.5 * (width2 + 1)))
        yTop    <- pmin(0.5, pmax(-0.5, y2 + 0.5 * (width2 - 1)))
        minY <- min(yBottom)
        maxY <- max(yTop)
        if (maxY == -0.5 || minY == 0.5) {
            ## outside drawing area
            return(grob(...)) # "empty" grob
        }
        x1 <- rep(c(-0.5, 0.5), each=2)
    }
    xList <- vector(mode = "list", length = nStripes)
    yList <- vector(mode = "list", length = nStripes)
    for (stripe in seq_len(nStripes)) {
        thisSin <- sinRot[stripe]
        thisCos <- cosRot[stripe]
        ## Rotation and scaling matrix: rotated full-size rectangle
        ## fills the whole viewport, clipping at the limit of the
        ## viewport (assuming default "npc" units, but other units are
        ## not really supported)
        R <- matrix(c(thisCos, thisSin, -thisSin, thisCos), 2, 2) *
            (abs(thisSin) + abs(thisCos))
        if (useX) {
            x1 <- c(xRight[stripe], xLeft[stripe],
                    xLeft[stripe], xRight[stripe])
            stripeVec <- drop(R %*% c(0, 1))
        } else {
            y1 <- c(yBottom[stripe], yTop[stripe],
                    yTop[stripe], yBottom[stripe])
            stripeVec <- drop(R %*% c(1, 0))
        }
        xy <- apply(R %*% rbind(x1, y1) + 0.5, 2, function(x) {
            if (isPerpendicular[stripe]) {
                mplier <- 0
                clipSide <- NA_real_ # 1=bottom, 2=left, 3=top, 4=right
            } else {
                clipSide <- which.max(pmax(0, c(-x[2], -x[1],
                                                x[2] - 1, x[1] - 1)))
                if (clipSide == 4) {
                    mplier <- (x[1] - 1) / stripeVec[1]
                } else if (clipSide == 2) {
                    mplier <- x[1] / stripeVec[1]
                } else if (clipSide == 3) {
                    mplier <- (x[2] - 1) / stripeVec[2]
                } else {
                    mplier <- x[2] / stripeVec[2]
                }
            }
            c(pmin(1, pmax(0, x - mplier * stripeVec)), clipSide)
        })
        lowPoints <- xy[1:2, 1:2]
        highPoints <- xy[1:2, 3:4]
        if (!isPerpendicular[stripe]) {
            clipSide <- xy[3, ]
            addLow <- clipSide[1] != clipSide[2]
            addHigh <- clipSide[3] != clipSide[4]
            if (addLow) {
                if (all(c(1, 2) %in% clipSide[1:2])) {
                    lowPoints <- cbind(lowPoints[, 1], c(0, 0), lowPoints[, 2])
                } else if (all(c(2, 3) %in% clipSide[1:2])) {
                    lowPoints <- cbind(lowPoints[, 1], c(0, 1), lowPoints[, 2])
                } else if (all(c(3, 4) %in% clipSide[1:2])) {
                    lowPoints <- cbind(lowPoints[, 1], c(1, 1), lowPoints[, 2])
                } else if (all(c(4, 1) %in% clipSide[1:2])) {
                    lowPoints <- cbind(lowPoints[, 1], c(1, 0), lowPoints[, 2])
                }
            }
            if (addHigh) {
                if (all(c(1, 2) %in% clipSide[3:4])) {
                    highPoints <-
                        cbind(highPoints[, 1], c(0, 0), highPoints[, 2])
                } else if (all(c(2, 3) %in% clipSide[3:4])) {
                    highPoints <-
                        cbind(highPoints[, 1], c(0, 1), highPoints[, 2])
                } else if (all(c(3, 4) %in% clipSide[3:4])) {
                    highPoints <-
                        cbind(highPoints[, 1], c(1, 1), highPoints[, 2])
                } else if (all(c(4, 1) %in% clipSide[3:4])) {
                    highPoints <-
                        cbind(highPoints[, 1], c(1, 0), highPoints[, 2])
                }
            }
            if (!(addLow || addHigh)) {
                if (identical(lowPoints[, 2], highPoints[, 1])) {
                    lowPoints <- lowPoints[, 1, drop=FALSE]
                } else if (identical(lowPoints[, 1],
                                     highPoints[, 2])) {
                    lowPoints <- lowPoints[, 2, drop=FALSE]
                }
            }
        }
        xList[[stripe]] <- c(lowPoints[1, ], highPoints[1, ]) *
            xScale2[stripe] + xShift2[stripe]
        yList[[stripe]] <- c(lowPoints[2, ], highPoints[2, ]) *
            yScale2[stripe] + yShift2[stripe]
    }
    polygonGrob(x = unlist(xList), y = unlist(yList),
                id.lengths = vapply(xList, length, 1), ...)
}
