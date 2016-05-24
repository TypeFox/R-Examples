### File R/dynTextGrob.R
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

## dynTextGrob is like a textGrob but with dynamic text size.  'name'
## and 'vp' are passed to grob(). 'label', 'x', 'y', 'just', 'hjust',
## 'vjust', 'rot' and 'gp' (sans "fontsize") are passed to textGrob().
##
## 'adjustJust' (a logical flag) toggles whether the justification
## settings passed to textGrob() will be adjusted (TRUE) or not
## (FALSE). The latter option can be useful for using previously
## batch-computed justification settings.
##
## If 'adjustJust' is TRUE, the [0, 1] range in the 'vjust' of
## plotmath labels will correspond to the area between the baseline
## (0) and baseline + one lineheight (1). Also, the 'vjust' of text
## labels will be adjusted to account for the dimensions of the
## plotmath labels (if any) and the (maximum) descent of the text
## labels themselves.  The adjustments, together with the fontsize
## computed for the labels, also ensure that the labels will fit in
## their designated space, for example when y == vjust and rot == 0.
##
## If the (recycled) element of 'rotJust' corresponding to a label is
## TRUE ('adjustJust' also required), justification is performed with
## respect to the bounding box of the rotated label.  If 'rotJust' is
## FALSE (the default), the label is justified traditionally, with
## respect to reading direction and the perpendicular direction, then
## rotated around the (x, y) point.
##
## 'rotVjust' and 'rotHjust' (NULL or recycled numeric vectors, with
## NA values allowed) can be used to adjust the justification of the
## labels with respect to each other when the corresponding element of
## 'rotJust' is TRUE. If NULL or NA, the justification settings used
## as the starting points before adjustment will be linearly
## interpolated according to rotation angle using the hjust and vjust
## settings that are arguably the most logical for the straight angles
## 0, 90, 180 and 270 degrees. For example, the default (NULL)
## behavior with a rotation of 90 degrees, hjust=1, vjust=1 and
## rotJust=TRUE is to align the bottom (last line) right corner of
## each label. With a rotation of 270 degrees but otherwise identical
## settings, the top (first line) left corners would be aligned. This
## dependence of aligment on rotation angle can be overridden with
## 'rotHjust' (reading direction) and 'rotVjust' (the perpendicular
## direction).
##
## If 'takeMeasurements' is TRUE, the function returns some
## measurements but does not return a grob. In this case, 'gp' can
## also be a list of "gpar" objects in list form.
##
## The number of labels 'n' is the maximum of the lengths of 'x' and
## 'y'.  The following vectors are (effectively) recycled to the same
## length: 'label', 'x', 'y', 'hjust', 'vjust', 'rot', 'rotJust'.
dynTextGrob <- function(label, x = 0.5, y = 0.5, width = 1, height = 1,
                        default.units = "npc", just = c(0.5, 0.5),
                        hjust = NULL, vjust = NULL, rot = 0, rotJust = TRUE,
                        rotHjust = NULL, rotVjust = NULL, resize = TRUE,
                        sizingWidth = NULL, sizingHeight = NULL,
                        adjustJust = TRUE, takeMeasurements = FALSE,
                        name = NULL, gp = gpar(), vp = NULL) {
    stopifnot(identical(resize, TRUE) || identical(resize, FALSE),
              length(length(label)) == 1,
              identical(adjustJust, TRUE) || identical(adjustJust, FALSE),
              identical(takeMeasurements, TRUE) ||
              identical(takeMeasurements, FALSE),
              is.null(sizingWidth) ||
              (is.numeric(sizingWidth) && length(sizingWidth) == 1 &&
               is.finite(sizingWidth) && sizingWidth >= 0),
              is.null(sizingHeight) ||
              (is.numeric(sizingHeight) && length(sizingHeight) == 1 &&
               is.finite(sizingHeight) && sizingHeight >= 0),
              is.logical(rotJust), length(rotJust) > 0, is.finite(rotJust),
              is.null(rotVjust) || is.numeric(rotVjust),
              is.null(rotHjust) || is.numeric(rotHjust),
              is.null(name) || (is.character(name) && length(name) == 1))
    REF_FONTSIZE <- 12
    if (is.null(sizingWidth) || is.null(sizingHeight) || !resize) {
        sizeGiven <- FALSE
    } else {
        sizeGiven <- TRUE
    }
    if (takeMeasurements) {
        n <- length(label)
        if (inherits(gp, "gpar")) {
            gpList <- gparToList(gp)
        } else {
            stopifnot(is.list(gp))
            gpList <- gp
        }
    } else {
        n <- max(length(x), length(y))
        stopifnot(inherits(gp, "gpar"))
        gpList <- gparToList(gp)
    }
    if (inherits(x, "unit")) {
        xUnit <- x
    } else {
        xUnit <- unit(x, default.units)
    }
    X <- rep(xUnit, length.out = n)
    if (inherits(y, "unit")) {
        yUnit <- y
    } else {
        yUnit <- unit(y, default.units)
    }
    Y <- rep(yUnit, length.out = n)
    xShift <- numeric(n)
    yShift <- numeric(n)
    if (inherits(width, "unit")) {
        width2 <- width
    } else {
        width2 <- unit(width, default.units)
    }
    if (inherits(height, "unit")) {
        height2 <- height
    } else {
        height2 <- unit(height, default.units)
    }
    tmp <- justHV(just, hjust, vjust, n)
    Hjust <- tmp[[1]]
    Vjust <- tmp[[2]]

    rot2 <- rep_len(rot %% 360, n)
    rotRad <- rot * (pi / 180)
    sinRot <- rep_len(sin(rotRad), n)
    cosRot <- rep_len(cos(rotRad), n)
    rotRad <- rep_len(rotRad, n)
    mul90 <- round(rot2 / 90) %% 4

    mfs <- plainLabelSize(label, n = n, rot = rot2, gpList = gpList,
                          resize = resize, refSize = REF_FONTSIZE)
    boxWidth <- mfs[["boxWidth"]]
    boxHeight <- mfs[["boxHeight"]]
    totalHeight <- mfs[["totalHeight"]]
    straightWidth <- mfs[["width"]]
    straightHeight <- mfs[["height"]]
    descent <- mfs[["descent"]]
    ascent <- mfs[["ascent"]]
    isCharacter <- mfs[["isCharacter"]]
    lSize <- mfs[["labelSize"]]
    nzDim <- lSize[["nzDim"]]

    if (adjustJust) {
        straightRot <- vapply(rot2,
                              function (x) {
                                  isTRUE(all.equal(45, abs(45 - x %% 90)))
                              }, TRUE)
        rotJust2 <- rep_len(rotJust, n)
        VjustOrig <- Vjust
        HjustOrig <- Hjust
        if (!is.null(rotVjust)) {
            rotVjust2 <- rep_len(rotVjust, n)
            rotVjustNA <- is.na(rotVjust2)
        } else {
            rotVjustNA <- rep.int(TRUE, n)
        }
        if (!is.null(rotHjust)) {
            rotHjust2 <- rep_len(rotHjust, n)
            rotHjustNA <- is.na(rotHjust2)
        } else {
            rotHjustNA <- rep.int(TRUE, n)
        }
        idxRotJust <- which(rotJust2)
        for (i in idxRotJust) {
            if (rotHjustNA[i]) {
                Hjust[i] <- approx(x = c(0, 90, 180, 270, 360),
                                   y = c(HjustOrig[i], VjustOrig[i],
                                       1 - HjustOrig[i], 1 - VjustOrig[i],
                                       HjustOrig[i]),
                                   xout = rot2[i], rule = 2)[["y"]]
            } else {
                Hjust[i] <- rotHjust2[i]
            }
            if (rotVjustNA[i]) {
                Vjust[i] <- approx(x = c(0, 90, 180, 270, 360),
                                   y = c(VjustOrig[i], 1 - HjustOrig[i],
                                       1 - VjustOrig[i], HjustOrig[i],
                                       VjustOrig[i]),
                                   xout = rot2[i], rule = 2)[["y"]]
            } else {
                Vjust[i] <- rotVjust2[i]
            }
        }

        jl <- justifyLabels(lSize, rot = rot2, hjust = Hjust, vjust = Vjust)
        bottomGap <- jl[["bottomGap"]]
        leftGap <- jl[["leftGap"]]
        bottomSpace <- rep.int(Inf, n)
        leftSpace <- rep.int(Inf, n)
        topSpace <- rep.int(Inf, n)
        rightSpace <- rep.int(Inf, n)
        maxWidth <- jl[["maxWidth"]]
        maxHeight <- jl[["maxHeight"]]
        paddedBoxWidth <- jl[["paddedBoxWidth"]]
        paddedBoxHeight <- jl[["paddedBoxHeight"]]

        nzRotJust <- idxRotJust[nzDim[idxRotJust]]
        nzStraightRot <- straightRot[nzRotJust]
        notRotJust <- !rotJust2
        Vjust[notRotJust] <- jl[["vjust"]][notRotJust]
        for (i in nzRotJust[!nzStraightRot]) {
            thisRot <- rot2[i]
            if (thisRot < 90) {
                bottomLeft <-
                    c((sinRot[i] * (maxHeight - bottomGap[i]) +
                       cosRot[i] * leftGap[i]) / paddedBoxWidth[i],
                      (cosRot[i] * bottomGap[i] +
                       sinRot[i] * leftGap[i]) / paddedBoxHeight[i])
                hUnit <-
                    c(cosRot[i] * straightWidth[i] / paddedBoxWidth[i],
                      sinRot[i] * straightWidth[i] / paddedBoxHeight[i])
                vTotal <- c(-sinRot[i] * totalHeight[i] / paddedBoxWidth[i],
                            cosRot[i] * totalHeight[i] / paddedBoxHeight[i])
                bottomSpace[i] <- bottomLeft[2] * paddedBoxHeight[i]
                leftSpace[i] <- (bottomLeft[1] + vTotal[1]) * paddedBoxWidth[i]
                topSpace[i] <- (1 - bottomLeft[2] - hUnit[2] - vTotal[2]) *
                    paddedBoxHeight[i]
                rightSpace[i] <- (1 - bottomLeft[1] - hUnit[1]) *
                    paddedBoxWidth[i]
            } else if (thisRot < 180) {
                sinTmp <- sin(rotRad[i] - 0.5 * pi)
                cosTmp <- sinRot[i]
                bottomLeft <-
                    c(1 - (sinTmp * leftGap[i] +
                           cosTmp * bottomGap[i]) / paddedBoxWidth[i],
                      (sinTmp * (maxHeight - bottomGap[i]) +
                       cosTmp * leftGap[i]) / paddedBoxHeight[i])
                hUnit <-
                    c(-sinTmp * straightWidth[i] / paddedBoxWidth[i],
                      cosTmp * straightWidth[i] / paddedBoxHeight[i])
                vTotal <- -c(cosTmp * totalHeight[i] / paddedBoxWidth[i],
                             sinTmp * totalHeight[i] / paddedBoxHeight[i])
                bottomSpace[i] <- (bottomLeft[2] + vTotal[2]) *
                    paddedBoxHeight[i]
                leftSpace[i] <- (bottomLeft[1] + hUnit[1] + vTotal[1]) *
                    paddedBoxWidth[i]
                topSpace[i] <- (1 - bottomLeft[2] - hUnit[2]) *
                    paddedBoxHeight[i]
                rightSpace[i] <- (1 - bottomLeft[1]) * paddedBoxWidth[i]
            } else if (thisRot < 270) {
                bottomLeft <-
                    c(1 + (sinRot[i] * (maxHeight - bottomGap[i]) +
                           cosRot[i] * leftGap[i]) / paddedBoxWidth[i],
                      1 + (cosRot[i] * bottomGap[i] +
                           sinRot[i] * leftGap[i]) / paddedBoxHeight[i])
                hUnit <-
                    c(cosRot[i] * straightWidth[i] / paddedBoxWidth[i],
                      sinRot[i] * straightWidth[i] / paddedBoxHeight[i])
                vTotal <- c(-sinRot[i] * totalHeight[i] / paddedBoxWidth[i],
                            cosRot[i] * totalHeight[i] / paddedBoxHeight[i])
                bottomSpace[i] <- (bottomLeft[2] + hUnit[2] + vTotal[2]) *
                    paddedBoxHeight[i]
                leftSpace[i] <- (bottomLeft[1] + hUnit[1]) * paddedBoxWidth[i]
                topSpace[i] <- (1 - bottomLeft[2]) * paddedBoxHeight[i]
                rightSpace[i] <- (1 - bottomLeft[1] - vTotal[1]) *
                    paddedBoxWidth[i]
            } else {
                sinTmp <- sin(rotRad[i] - 1.5 * pi)
                cosTmp <- -sinRot[i]
                bottomLeft <-
                    c((sinTmp * leftGap[i] +
                       cosTmp * bottomGap[i]) / paddedBoxWidth[i],
                      1 - (sinTmp * (maxHeight - bottomGap[i]) +
                           cosTmp * leftGap[i]) / paddedBoxHeight[i])
                hUnit <-
                    c(sinTmp * straightWidth[i] / paddedBoxWidth[i],
                      -cosTmp * straightWidth[i] / paddedBoxHeight[i])
                vTotal <- c(cosTmp * totalHeight[i] / paddedBoxWidth[i],
                            sinTmp * totalHeight[i] / paddedBoxHeight[i])
                bottomSpace[i] <- (bottomLeft[2] + hUnit[2]) *
                    paddedBoxHeight[i]
                leftSpace[i] <- bottomLeft[1] * paddedBoxWidth[i]
                topSpace[i] <- (1 - bottomLeft[2] - vTotal[2]) *
                    paddedBoxHeight[i]
                rightSpace[i] <- (1 - bottomLeft[1] - hUnit[1] - vTotal[1]) *
                    paddedBoxWidth[i]
            }
            if (isCharacter[i]) {
                vUnit <- vTotal * (straightHeight[i] / totalHeight[i])
                origin <- bottomLeft +
                    (descent[i] / totalHeight[i]) * vTotal
            } else {
                vUnit <- vTotal
                origin <- bottomLeft
            }
            A <- hUnit[1]
            B <- vUnit[1]
            C <- hUnit[2]
            D <- vUnit[2]
            a <- HjustOrig[i] - origin[1]
            b <- VjustOrig[i] - origin[2]
            hj <- (b * B - a * D) / (C * B - A * D)
            Vjust[i] <- (b * A - a * C) / (D * A - B * C)
            ## It's necessary to adjust the X, Y coordinates in order
            ## to achieve a horizontal shift. Adjusting Hjust directly
            ## would throw off the justification of multi-line text
            ## labels.
            hShift <- (Hjust[i] - hj) * straightWidth[i]
            xShift[i] <- cosRot[i] * hShift
            yShift[i] <- sinRot[i] * hShift
        }
        for (i in nzRotJust[nzStraightRot]) {
            thisMul <- mul90[i]
            if (thisMul == 0) {
                bottomLeft <- c(leftGap[i] / maxWidth,
                                bottomGap[i] / maxHeight)
                hUnit <- c(straightWidth[i] / maxWidth, 0)
                vTotal <- c(0, totalHeight[i] / maxHeight)
                bottomSpace[i] <- bottomGap[i]
                leftSpace[i] <- leftGap[1]
                topSpace[i] <- maxHeight - bottomGap[i] - totalHeight[i]
                rightSpace[i] <- maxWidth - leftGap[i] - straightWidth[i]
            } else if (thisMul == 1) {
                bottomLeft <- c(1 - bottomGap[i] / maxHeight,
                                leftGap[i] / maxWidth)
                hUnit <- c(0, straightWidth[i] / maxWidth)
                vTotal <- c(-totalHeight[i] / maxHeight, 0)
                bottomSpace[i] <- leftGap[i]
                leftSpace[i] <- maxHeight - bottomGap[i] - totalHeight[i]
                topSpace[i] <- maxWidth - leftGap[i] - straightWidth[i]
                rightSpace[i] <- bottomGap[i]
            } else if (thisMul == 2) {
                bottomLeft <- 1 - c(leftGap[i] / maxWidth,
                                    bottomGap[i] / maxHeight)
                hUnit <- c(-straightWidth[i] / maxWidth, 0)
                vTotal <- c(0, -totalHeight[i] / maxHeight)
                bottomSpace[i] <- maxHeight - bottomGap[i] - totalHeight[i]
                leftSpace[i] <- maxWidth - leftGap[i] - straightWidth[i]
                topSpace[i] <- bottomGap[i]
                rightSpace[i] <- leftGap[i]
            } else {
                bottomLeft <- c(bottomGap[i] / maxHeight,
                                1 - leftGap[i] / maxWidth)
                hUnit <- c(0, -straightWidth[i] / maxWidth)
                vTotal <- c(totalHeight[i] / maxHeight, 0)
                bottomSpace[i] <- maxWidth - leftGap[i] - straightWidth[i]
                leftSpace[i] <- bottomGap[i]
                topSpace[i] <- leftGap[i]
                rightSpace[i] <- maxHeight - bottomGap[i] - totalHeight[i]
            }
            if (isCharacter[i]) {
                vUnit <- vTotal * (straightHeight[i] / totalHeight[i])
                origin <- bottomLeft +
                    (descent[i] / totalHeight[i]) * vTotal
            } else {
                vUnit <- vTotal
                origin <- bottomLeft
            }
            if (thisMul == 0 || thisMul == 2) {
                Vjust[i] <- (VjustOrig[i] - origin[2]) / vUnit[2]
                hj <- (HjustOrig[i] - origin[1]) / hUnit[1]
            } else {
                Vjust[i] <- (HjustOrig[i] - origin[1]) / vUnit[1]
                hj <- (VjustOrig[i] - origin[2]) / hUnit[2]
            }
            hShift <- (Hjust[i] - hj) * straightWidth[i]
            xShift[i] <- cosRot[i] * hShift
            yShift[i] <- sinRot[i] * hShift
        }

        bottomSpace <- min(bottomSpace)
        leftSpace <- min(leftSpace)
        topSpace <- min(topSpace)
        rightSpace <- min(rightSpace)
        if (!is.finite(bottomSpace)) {
            bottomSpace <- 0
        }
        if (!is.finite(leftSpace)) {
            leftSpace <- 0
        }
        if (!is.finite(topSpace)) {
            topSpace <- 0
        }
        if (!is.finite(rightSpace)) {
            rightSpace <- 0
        }
        xShift <- xShift + HjustOrig * (leftSpace + rightSpace) - leftSpace
        yShift <- yShift + VjustOrig * (bottomSpace + topSpace) - bottomSpace

        if (sizeGiven) {
            sizingWidth2 <- sizingWidth
            sizingHeight2 <- sizingHeight
        } else {
            sizingWidth2 <- max(paddedBoxWidth) - leftSpace - rightSpace
            sizingHeight2 <- max(paddedBoxHeight) - bottomSpace - topSpace
        }
    } else if (sizeGiven) {
        sizingWidth2 <- sizingWidth
        sizingHeight2 <- sizingHeight
    } else {
        sizingWidth2 <- max(boxWidth)
        sizingHeight2 <- max(boxHeight)
    }

    if (takeMeasurements) {
        if (!resize) {
            leftRightBottomTop <-
                border2.dynText(Hjust, Vjust, sinRot, cosRot,
                                straightWidth, straightHeight,
                                ascent, descent, isCharacter)
            leftCoords <- leftRightBottomTop[["left"]] + xShift
            rightCoords <- leftRightBottomTop[["right"]] + xShift
            bottomCoords <- leftRightBottomTop[["bottom"]] + yShift
            topCoords <- leftRightBottomTop[["top"]] + yShift
            return(list(sizingWidth = sizingWidth2, isCharacter = isCharacter,
                        sizingHeight = sizingHeight2, vjust = Vjust,
                        hjust = Hjust, nzDim = nzDim, left = leftCoords,
                        right = rightCoords, bottom = bottomCoords,
                        top = topCoords, xShift = xShift, yShift = yShift,
                        straightWidth = straightWidth, ascent = ascent,
                        straightHeight = straightHeight, descent = descent))
        } else {
            return(list(sizingWidth = sizingWidth2, isCharacter = isCharacter,
                        sizingHeight = sizingHeight2, nzDim = nzDim))
        }
    }


    if (resize) {
        gpList["fontsize"] <- NULL
        gp2 <- do.call("gpar", gpList)
        gpList[["cex"]] <- 1
        X2 <- X
        Y2 <- Y
        vpCex <- NULL
    } else {
        ## "Not resize" is implemented by setting missing gp settings
        ## to current values.  Additionally, cex needs to be recorded
        ## and checked at drawing time (makeContent.dynText).
        currentGP <- gparToList(get.gpar())
        vpCex <- currentGP[["cex"]]
        if (is.null(vpCex)) {
            vpCex <- 1
        } else {
            currentGP[["cex"]] <- 1
        }
        currentGP[names(gpList)] <- gpList
        gpList <- gparToList(currentGP)
        gp2 <- do.call("gpar", gpList)
        X2 <- X + unit(xShift, "inches")
        Y2 <- Y + unit(yShift, "inches")
    }
    getFontsize <-
        makeGetFontsize.dynText(label, gpList, sizingWidth2, sizingHeight2,
                                REF_FONTSIZE, n, width2, height2, resize)
    tmp <- setLabel.dynText(label, n, isCharacter, X2, Y2, Hjust, Vjust, rot2)
    tgList <- tmp[[1]]
    charIdx <- tmp[[2]]
    otherIdx <- tmp[[3]]
    oneGrob <- tmp[[4]]

    gTree(children = tgList,
          internals = list(x = X, y = Y, xShift = xShift, yShift = yShift,
              oneGrob = oneGrob, sizingWidth = sizingWidth2, vpCex = vpCex,
              sizingHeight = sizingHeight2, widthOrig = sizingWidth2,
              heightOrig = sizingHeight2, charIdx = charIdx,
              otherIdx = otherIdx, hjust = Hjust, vjust = Vjust, rot = rot2,
              sinRot = sinRot, cosRot = cosRot, ascent = ascent,
              descent = descent, straightWidth = straightWidth,
              straightHeight = straightHeight, refSize = REF_FONTSIZE,
              getFontsize = getFontsize),
          label = label, x = xUnit, y = yUnit,
          width = width2, height = height2, just = just, hjust = hjust,
          vjust = vjust, sizingWidth = sizingWidth, sizingHeight = sizingHeight,
          rot = rot, rotJust = rotJust, rotVjust = rotVjust,
          rotHjust = rotHjust, adjustJust = adjustJust, resize = resize,
          name = name, gp = gp2, vp = vp, cl = "dynText")
}

### Auxiliary functions for dynText

labelSize <- function(label, gpList, fontsize, resize) {
    n <- length(label)
    lineHeight <- numeric(n)
    width <- numeric(n)
    height <- numeric(n)
    ascent <- numeric(n)
    descent <- numeric(n)
    labelIsList <- is.list(label)
    for (i in seq_len(n)) {
        pushViewport(viewport(gp = do.call("gpar", gpList[[i]])),
                     recording = FALSE)
        lineHeight[i] <- convertUnit(stringHeight("x"),
                                     axisFrom = "y", valueOnly = TRUE,
                                     unitTo = "inches", typeFrom = "dimension")
        if (labelIsList) {
            tmpGrob <- textGrob(label[[i]])
        } else {
            tmpGrob <- textGrob(label[i])
        }
        width[i] <- convertUnit(grobWidth(tmpGrob),
                                axisFrom = "x", valueOnly = TRUE,
                                unitTo = "inches", typeFrom = "dimension")
        height[i] <- convertUnit(grobHeight(tmpGrob),
                                 axisFrom = "y", valueOnly = TRUE,
                                 unitTo = "inches", typeFrom = "dimension")
        ascent[i] <- convertUnit(grobAscent(tmpGrob),
                                 axisFrom = "y", valueOnly = TRUE,
                                 unitTo = "inches", typeFrom = "dimension")
        descent[i] <- convertUnit(grobDescent(tmpGrob),
                                  axisFrom = "y", valueOnly = TRUE,
                                  unitTo = "inches", typeFrom = "dimension")
        popViewport(recording = FALSE)
    }
    isCharacter <- if (is.language(label)) {
        rep.int(FALSE, length(label))
    } else if (is.list(label)) {
        !vapply(label, is.language, TRUE)
    } else {
        rep.int(TRUE, length(label))
    }
    if (resize) {
        width <- width / fontsize
        height <- height / fontsize
        ascent <- ascent / fontsize
        descent <- descent / fontsize
        lineHeight <- lineHeight / fontsize
    }

    nzDim <- width > 0 & height > 0
    nzIdx <- which(nzDim)
    list(width = width, height = height, ascent = ascent, descent = descent,
         isCharacter = isCharacter, lineHeight = lineHeight, nzDim = nzDim,
         nzIdx = nzIdx)
}

## Maximum width and height of 'label'.  Returns numeric values:
## width and height (points) per one fontsize unit (points).  The idea
## is to test if the proportions of a textgrob are linear with respect
## to fontsize.  A small selection of different font sizes is used for
## testing.  If resize is FALSE, computes actual fixed sizes in inches.
plainLabelSize <- function(label, rot = 0, n = length(label),
                           gpList = list(list()), resize = TRUE, refSize) {
    if (length(gpList) > 0 && all(vapply(gpList, is.list, TRUE))) {
        gpList2 <- rep_len(gpList, n)
    } else {
        gpList2 <- rep_len(list(gpList), n)
    }
    repLabel <- rep(label, length.out = n)
    if (resize) {
        ## Test all labels with a reference font size
        fontsize <- refSize
        for (i in seq_len(n)) {
            ## Assume cex 1 => allow different cex to have an effect later
            gpList2[[i]][["cex"]] <- 1
            gpList2[[i]][["fontsize"]] <- fontsize
        }
    } else {
        fontsize <- gpList2[[1]][["fontsize"]]
    }
    lSize <- labelSize(repLabel, gpList2, fontsize, resize)
    width <- lSize[["width"]]
    height <- lSize[["height"]]
    ascent <- lSize[["ascent"]]
    descent <- lSize[["descent"]]
    isCharacter <- lSize[["isCharacter"]]
    nzIdx <- lSize[["nzIdx"]]

    rot2 <- rep_len(rot * (pi / 180), n)[nzIdx]
    absSinRot2 <- abs(sin(rot2))
    absCosRot2 <- abs(cos(rot2))
    totalHeight <- height
    totalHeight[isCharacter] <- totalHeight[isCharacter] + descent[isCharacter]
    boxWidth <- numeric(n)
    boxHeight <- numeric(n)
    tmpWidth <- width[nzIdx]
    tmpHeight <- totalHeight[nzIdx]
    boxWidth1 <- absSinRot2 * tmpHeight + absCosRot2 * tmpWidth
    boxHeight1 <- absSinRot2 * tmpWidth + absCosRot2 * tmpHeight
    boxWidth[nzIdx] <- boxWidth1
    boxHeight[nzIdx] <- boxHeight1

    list(boxWidth = boxWidth, boxHeight = boxHeight,
         width = width, height = height, totalHeight = totalHeight,
         ascent = ascent, descent = descent, fontsize = fontsize,
         isCharacter = isCharacter, labelSize = lSize)
}

justifyLabels <- function(lSize, rot = 0, vjust = 0.5, hjust = 0.5) {
    width <- lSize[["width"]]
    height <- lSize[["height"]]
    ascent <- lSize[["ascent"]]
    descent <- lSize[["descent"]]
    lineHeight <- lSize[["lineHeight"]]
    isCharacter <- lSize[["isCharacter"]]
    nzIdx <- lSize[["nzIdx"]]
    n <- length(width)
    nzIsCharacter <- isCharacter[nzIdx]
    idxPlotmath <- nzIdx[!nzIsCharacter]
    idxCharacter <- nzIdx[nzIsCharacter]
    vjust2 <- rep_len(vjust, n)
    hjust2 <- rep_len(hjust, n)

    maxDescent <- max(descent)
    aPlotmath <- ascent[idxPlotmath]
    dPlotmath <- descent[idxPlotmath]
    hPlotmath <- height[idxPlotmath]
    lhPlotmath <- lineHeight[idxPlotmath]
    hChar <- height[idxCharacter]
    topHang <- max(0, aPlotmath - lhPlotmath)
    vjustPlotmath <-
        (dPlotmath - maxDescent + vjust2[idxPlotmath] *
         (max(0, aPlotmath, lhPlotmath) + maxDescent)) / hPlotmath
    vjust2[idxPlotmath] <- vjustPlotmath
    vjustChar <- (-maxDescent + vjust2[idxCharacter] *
                  (maxDescent + topHang + hChar)) / hChar
    vjust2[idxCharacter] <- vjustChar
    charTop <- (1 - vjustChar) * hChar
    plotmathTop <- (1 - vjustPlotmath) * hPlotmath
    charBottom <- -vjustChar * hChar - descent[idxCharacter]
    plotmathBottom <- -vjustPlotmath * hPlotmath
    maxTop <- max(plotmathTop, charTop, -Inf)
    minBottom <- min(plotmathBottom, charBottom, Inf)
    bottomCoord <- numeric(n)
    topCoord <- numeric(n)
    bottomCoord[idxCharacter] <- charBottom
    topCoord[idxCharacter] <- charTop
    bottomCoord[idxPlotmath] <- plotmathBottom
    topCoord[idxPlotmath] <- plotmathTop

    maxWidth <- max(width)
    maxHeight <- maxTop - minBottom
    if (!is.finite(maxHeight)) {
        maxHeight <- 0
    }
    leftGap <- numeric(n)
    leftGap[nzIdx] <- hjust2[nzIdx] * (maxWidth - width[nzIdx])
    rot2 <- rep_len(rot * (pi / 180), n)[nzIdx]
    absSinRot2 <- abs(sin(rot2))
    absCosRot2 <- abs(cos(rot2))
    bottomGap <- numeric(n)
    bottomGap[nzIdx] <- bottomCoord[nzIdx] - minBottom

    paddedBoxWidth <- numeric(n)
    paddedBoxHeight <- numeric(n)
    paddedBoxWidth[nzIdx] <- absSinRot2 * maxHeight + absCosRot2 * maxWidth
    paddedBoxHeight[nzIdx] <- absSinRot2 * maxWidth + absCosRot2 * maxHeight

    list(paddedBoxWidth = paddedBoxWidth, paddedBoxHeight = paddedBoxHeight,
         vjust = vjust2, maxHeight = maxHeight, maxWidth = maxWidth,
         bottomGap = bottomGap, leftGap = leftGap)
}

## Construct one or two "text" grobs holding the contents
setLabel.dynText <- function(label, n, isCharacter, x, y, hjust, vjust, rot) {
    if (is.list(label)) {
        nCharacter <- sum(isCharacter)
        labelRep <- rep(label, length.out = n)
        if (nCharacter == n) {
            labelChar <- vapply(labelRep, function(x) x[[1]], character(1))
            tg <- textGrob(label = labelChar, x = x, y = y,
                           hjust = hjust, vjust = vjust, rot = rot,
                           name = "one")
            list(gList(tg), 1:n, integer(0), TRUE)
        } else if (nCharacter == 0) {
            labelOther <- do.call("c", labelRep)
            tg <- textGrob(label = labelOther, x = x, y = y,
                           hjust = hjust, vjust = vjust, rot = rot,
                           name = "one")
            list(gList(tg), integer(0), 1:n, TRUE)
        } else {
            charIdx <- which(isCharacter)
            labelChar <- vapply(labelRep[charIdx],
                                function(x) x[[1]], character(1))
            tgChar <-
                textGrob(label = labelChar,
                         x = x[charIdx], y = y[charIdx],
                         hjust = hjust[charIdx], vjust = vjust[charIdx],
                         rot = rot[charIdx], name = "char")
            otherIdx <- which(!isCharacter)
            labelOther <- do.call("c", labelRep[otherIdx])
            tgOther <-
                textGrob(label = labelOther,
                         x = x[otherIdx], y = y[otherIdx],
                         hjust = hjust[otherIdx], vjust = vjust[otherIdx],
                         rot = rot[otherIdx], name = "other")
            list(gList(tgChar, tgOther), charIdx, otherIdx, FALSE)
        }
    } else {
        tg <- textGrob(label = label, x = x, y = y,
                       hjust = hjust, vjust = vjust, rot = rot,
                       name = "one")
        if (is.language(label)) {
            list(gList(tg), integer(0), 1:n, TRUE)
        } else {
            list(gList(tg), 1:n, integer(0), TRUE)
        }
    }
}

border.dynText <- function(x, chull = FALSE) {
    internals <- makeContent.dynText(x)[["internals"]]
    X <- internals[["x"]] + unit(internals[["xShift"]], "inches")
    Y <- internals[["y"]] + unit(internals[["yShift"]], "inches")
    X <- convertUnit(X, unitTo = "inches", axisFrom = "x",
                     typeFrom = "location", valueOnly = TRUE)
    Y <- convertUnit(Y, unitTo = "inches", axisFrom = "y",
                     typeFrom = "location", valueOnly = TRUE)
    hjust <- internals[["hjust"]]
    vjust <- internals[["vjust"]]
    sinRot <- internals[["sinRot"]]
    cosRot <- internals[["cosRot"]]
    straightWidth <- internals[["straightWidth"]]
    straightHeight <- internals[["straightHeight"]]
    width <- straightWidth
    height <- straightHeight
    if (chull) {
        ascent <- internals[["ascent"]]
        descent <- internals[["descent"]]
        otherIdx <- internals[["otherIdx"]]
        ascent[otherIdx] <- height[otherIdx]
        descent[otherIdx] <- 0
    }

    tmp <- width * cosRot
    xLeft <- X - tmp * hjust
    xRight <- xLeft + tmp
    xBottom <- height * vjust * sinRot
    if (chull) {
        xTop <- xBottom - ascent * sinRot
        xBottom <- xBottom + descent * sinRot
    } else {
        xTop <- xBottom - height * sinRot
    }
    xTopLeft <- xTop + xLeft
    xBottomLeft <- xBottom + xLeft
    xTopRight <- xTop + xRight
    xBottomRight <- xBottom + xRight

    tmp <- width * sinRot
    yLeft <- Y - tmp * hjust
    yRight <- yLeft + tmp
    yBottom <- -height * vjust * cosRot
    if (chull) {
        yTop <- yBottom + ascent * cosRot
        yBottom <- yBottom - descent * cosRot
    } else {
        yTop <- yBottom + height * cosRot
    }
    yTopLeft <- yTop + yLeft
    yBottomLeft <- yBottom + yLeft
    yTopRight <- yTop + yRight
    yBottomRight <- yBottom + yRight

    if (chull) {
        ## Convex hull: x coords in 1st column, y coords in 2nd column
        xy <- cbind(c(xTopLeft, xTopRight, xBottomLeft, xBottomRight),
                    c(yTopLeft, yTopRight, yBottomLeft, yBottomRight))
        xy[chull(xy), ]
    } else {
        xRange <- range(xTopLeft, xTopRight, xBottomLeft, xBottomRight)
        yRange <- range(yTopLeft, yTopRight, yBottomLeft, yBottomRight)
        ## Coords of bounding box: c(left, right, bottom, top)
        c(xRange, yRange)
    }
}

border2.dynText <- function(hjust, vjust, sinRot, cosRot,
                            width, height, ascent, descent,
                            isCharacter) {
    ascent2 <- ascent
    descent2 <- descent
    isNotCharacter <- which(!isCharacter)
    ascent2[isNotCharacter] <- height[isNotCharacter]
    descent2[isNotCharacter] <- 0

    tmp <- width * cosRot
    xLeft <- -tmp * hjust
    xRight <- xLeft + tmp
    tmp <- height * vjust * sinRot
    xTop <- tmp - ascent2 * sinRot
    xBottom <- tmp + descent2 * sinRot
    xTopLeft <- xTop + xLeft
    xBottomLeft <- xBottom + xLeft
    xTopRight <- xTop + xRight
    xBottomRight <- xBottom + xRight

    tmp <- width * sinRot
    yLeft <- -tmp * hjust
    yRight <- yLeft + tmp
    tmp <- -height * vjust * cosRot
    yTop <- tmp + ascent2 * cosRot
    yBottom <- tmp - descent2 * cosRot
    yTopLeft <- yTop + yLeft
    yBottomLeft <- yBottom + yLeft
    yTopRight <- yTop + yRight
    yBottomRight <- yBottom + yRight

    list(left = pmin(xTopLeft, xTopRight, xBottomLeft, xBottomRight),
         right = pmax(xTopLeft, xTopRight, xBottomLeft, xBottomRight),
         bottom = pmin(yTopLeft, yTopRight, yBottomLeft, yBottomRight),
         top = pmax(yTopLeft, yTopRight, yBottomLeft, yBottomRight))
}

makeGetFontsize.dynText <- function(label, gpList, sizingWidth, sizingHeight,
                                    refSize, n, width, height, resize) {
    if (resize && sizingWidth > 0 && sizingHeight > 0) {
        ## Test one label with a selection of arbitrary font sizes.
        ## If there are several labels, "Dummy" is used. Otherwise,
        ## the single unique label is tested.
        FS_MIN_LOG2 <- 1  # min font size 2
        FS_MAX_LOG2 <- 10 # max font size 1024
        N_SEQ <- 37
        fsSeq <- sort(union(2^seq(from = FS_MIN_LOG2, to = FS_MAX_LOG2,
                                  length.out = N_SEQ), refSize))
        nSizes <- length(fsSeq) # N_SEQ or N_SEQ + 1
        nLabel <- length(label)
        if (n == 1 || nLabel == 1 ||
            (!is.language(label) && length(unique(label)) == 1)) {
            if (nLabel > 1) {
                testLabel <- label[[1]]
            } else {
                testLabel <- label
            }
            if (length(testLabel) > 1) {
                testLabel <- testLabel[1]
            }
        } else {
            testLabel <- "Dummy"
        }
        idxRef <- which(fsSeq == refSize)[1]
        testWidth <- numeric(nSizes)
        testHeight <- numeric(nSizes)
        gpList2 <- gpList
        for (i in seq_len(nSizes)) {
            fontsize <- fsSeq[[i]]
            gpList2[["fontsize"]] <- fontsize
            tmp <- labelSize(testLabel, list(gpList2), fontsize, TRUE)
            testWidth[i] <- tmp[["width"]]
            testHeight[i] <- tmp[["ascent"]] + tmp[["descent"]]
        }
        if (all(c(testWidth, testHeight) > 0)) {
            scaleTableX <- testWidth * (sizingWidth / testWidth[idxRef])
            scaleTableY <- testHeight * (sizingHeight / testHeight[idxRef])
        } else {
            scaleTableX <- fsSeq * (sizingWidth / refSize)
            scaleTableY <- fsSeq * (sizingHeight / refSize)
        }
        scaleTableX <- scaleTableX + sd(scaleTableX)
        scaleTableY <- scaleTableY + sd(scaleTableY)
        ## xApprox and yApprox are created in a complicated manner in
        ## order to avoid spurious checkUsagePackage(all=TRUE)
        ## warnings
        xApprox <- function(v) { v }
        yApprox <- function(v) { v }
        xApproxTmp <-
            approxfun(x = fsSeq, y = scaleTableX, method = "linear", rule = 2)
        yApproxTmp <-
            approxfun(x = fsSeq, y = scaleTableY, method = "linear", rule = 2)
        formals(xApprox) <- formals(xApproxTmp)
        formals(yApprox) <- formals(yApproxTmp)
        body(xApprox) <- body(xApproxTmp)
        body(yApprox) <- body(yApproxTmp)
        environment(xApprox) <- environment(xApproxTmp)
        environment(yApprox) <- environment(yApproxTmp)
        ## Functions to minimize when optimizing fontsize with respect
        ## to space available in the x and y directions
        optimX <- function(fontsize, widthInches) {
            (xApprox(fontsize) * fontsize - widthInches)^2
        }
        optimY <- function(fontsize, heightInches) {
            (yApprox(fontsize) * fontsize - heightInches)^2
        }
        getFontsize <- function(...) {
            widthInches <-
                convertUnit(width, axisFrom = "x",
                            unitTo="inches", typeFrom = "dimension",
                            valueOnly=TRUE)
            heightInches <-
                convertUnit(height, axisFrom = "y",
                            unitTo="inches", typeFrom = "dimension",
                            valueOnly=TRUE)
            xFS <- widthInches / sizingWidth
            yFS <- heightInches / sizingHeight
            min(optimize(optimX, widthInches = widthInches,
                         maximum = FALSE, lower = xFS * 0.8,
                         upper = xFS * 1.2)[["minimum"]],
                optimize(optimY, heightInches = heightInches,
                         maximum = FALSE, lower = yFS * 0.8,
                         upper = yFS * 1.2)[["minimum"]])
        }
        ## Making the environment of getFontsize() as small as possible
        newEnv <- new.env(hash = FALSE, parent = environment(sisal))
        assign("width", width, pos = newEnv)
        assign("height", height, pos = newEnv)
        assign("sizingWidth", sizingWidth, pos = newEnv)
        assign("sizingHeight", sizingHeight, pos = newEnv)
        assign("xApprox", xApprox, pos = newEnv)
        assign("yApprox", yApprox, pos = newEnv)
        environment(optimX) <- newEnv
        assign("optimX", optimX, pos = newEnv)
        environment(optimY) <- newEnv
        assign("optimY", optimY, pos = newEnv)
        environment(getFontsize) <- newEnv
        getFontsize
    } else {
        fontsize <- gpList[["fontsize"]]
        if (is.null(fontsize)) {
            getFontsize <- function(...) {
                unname(unlist(get.gpar("fontsize")))
            }
            environment(getFontsize) <- environment(sisal)
            getFontsize
        } else {
            getFontsize <- function(...) {}
            body(getFontsize) <- fontsize
            environment(getFontsize) <- emptyenv()
            getFontsize
        }
    }
}

### S3 methods for dynText

widthDetails.dynText <- function(x, ...) {
    xy <- border.dynText(x, chull = FALSE)
    unit(xy[2] - xy[1], "inches")
}

heightDetails.dynText <- function(x, ...) {
    xy <- border.dynText(x, chull = FALSE)
    unit(xy[4] - xy[3], "inches")
}

ascentDetails.dynText <- function(x, ...) {
    label <- x[["label"]]
    resize <- x[["resize"]]
    if (length(label) == 1) {
        if (resize) {
            gp <- gpar(fontsize = x[["internals"]][["getFontsize"]]())
            pushViewport(viewport(gp = gp), recording = FALSE)
        }
        if (is.list(label)) {
            ascent <- calcStringMetric(label[[1]])[["ascent"]]
        } else {
            ascent <- calcStringMetric(label)[["ascent"]]
        }
        if (resize) {
            popViewport(recording = FALSE)
        }
        unit(ascent, "inches")
    } else {
        heightDetails.dynText(x)
    }
}

descentDetails.dynText <- function(x, ...) {
    label <- x[["label"]]
    resize <- x[["resize"]]
    if (length(label) == 1) {
        if (resize) {
            gp <- gpar(fontsize = x[["internals"]][["getFontsize"]]())
            pushViewport(viewport(gp = gp), recording = FALSE)
        }
        if (is.list(label)) {
            descent <- calcStringMetric(label[[1]])[["descent"]]
        } else {
            descent <- calcStringMetric(label)[["descent"]]
        }
        if (resize) {
            popViewport(recording = FALSE)
        }
        unit(descent, "inches")
    } else {
        unit(0, "inches")
    }
}

## xDetails and yDetails methods of dynText are vectorized
## w.r.t. theta. It would be nice if grid could take advantage of this
## because repeated calls to border.dynText() can be very slow.
xDetails.dynText <- function(x, theta, ...) {
    xy <- border.dynText(x, chull = TRUE)
    pGrob <- polygonGrob(x = xy[, 1], y = xy[, 2], default.units = "inches")
    do.call("unit.c", lapply(theta, function (x) xDetails(pGrob, x)))
}

yDetails.dynText <- function(x, theta, ...) {
    xy <- border.dynText(x, chull = TRUE)
    pGrob <- polygonGrob(x = xy[, 1], y = xy[, 2], default.units = "inches")
    do.call("unit.c", lapply(theta, function (x) yDetails(pGrob, x)))
}

makeContent.dynText <- function(x, ...) {
    if (x[["resize"]]) {
        internals <- x[["internals"]]
        fontsize <- internals[["getFontsize"]]()
        gp <- gpar(fontsize = fontsize)
        x2 <- x
        if (x[["adjustJust"]]) {
            gpList <- gparToList(x[["gp"]])
            gpList[["fontsize"]] <- fontsize
            measurements <-
                dynTextGrob(label = x[["label"]], just = x[["just"]],
                            hjust = x[["hjust"]], vjust = x[["vjust"]],
                            rot = x[["rot"]], rotJust = x[["rotJust"]],
                            rotHjust = x[["rotHjust"]],
                            rotVjust = x[["rotVjust"]],
                            resize = FALSE, adjustJust = TRUE,
                            takeMeasurements = TRUE,
                            gp = do.call("gpar", gpList))
            slotNames <- c("xShift", "yShift", "hjust", "vjust", "ascent",
                           "descent", "sizingWidth", "sizingHeight",
                           "straightWidth", "straightHeight")
            internals[slotNames]<- measurements[slotNames]
            X <- internals[["x"]] + unit(measurements[["xShift"]], "inches")
            Y <- internals[["y"]] + unit(measurements[["yShift"]], "inches")
            hjust <- measurements[["hjust"]]
            vjust <- measurements[["vjust"]]
            f <- function(...) {}
            body(f) <- fontsize
            environment(f) <- emptyenv()
            internals[["getFontsize"]] <- f
            x2[["internals"]] <- internals
            if (internals[["oneGrob"]]) {
                editGrob(x2, "one", x = X, y = Y, hjust = hjust, vjust = vjust,
                         gp = gp)
            } else {
                charIdx <- internals[["charIdx"]]
                otherIdx <- internals[["otherIdx"]]
                x2 <- editGrob(x2, "char", x = X[charIdx], y = Y[charIdx],
                               hjust = hjust[charIdx], vjust = vjust[charIdx],
                               gp = gp)
                editGrob(x2, "other", x = X[otherIdx], y = Y[otherIdx],
                         hjust = hjust[otherIdx], vjust = vjust[otherIdx],
                         gp = gp)
            }
        } else {
            for (child in childNames(x)) {
                x2 <- editGrob(x2, child, gp = gp)
            }
            x2
        }
    } else {
        currentCex <- tryCatch(unname(unlist(get.gpar("cex"))),
                               error = function (...) 1)
        vpCex <- x[["internals"]][["vpCex"]]
        if (isTRUE(all.equal(vpCex, currentCex))) {
            x
        } else {
            gp <- gpar(cex = vpCex / currentCex)
            x2 <- x
            for (child in childNames(x)) {
                x2 <- editGrob(x2, child, gp = gp)
            }
            x2
        }
    }
}

validDetails.dynText <- function(x, ...) {
    slotNames <- names(x)
    requiredNames <-
        c("internals", "label", "x", "y", "width",
          "height", "just", "hjust", "vjust", "sizingWidth", "sizingHeight",
          "rot", "rotJust", "rotVjust", "rotHjust", "adjustJust", "resize",
          "name", "gp", "vp")

    ## Testing that all slots are present
    if (!all(requiredNames %in% slotNames)) {
        stop("slots missing")
    }

    ## Testing 'resize'
    resize <- x[["resize"]]
    if (!(identical(resize, TRUE) || identical(resize, FALSE))) {
        stop(gettextf("'%s' must be a logical flag", "resize",
                      domain = "R-sisal"), domain = NA)
    }

    ## Testing internal variables
    internals <- x[["internals"]]
    if (!is.list(internals)) {
        stop(gettextf("'%s' must be a list", 'x[["internals"]]',
                      domain = "R-sisal"), domain = NA)
    }
    varNames <- c("x", "y", "oneGrob", "sizingWidth", "sizingHeight",
                  "widthOrig", "heightOrig", "charIdx", "otherIdx",
                  "hjust", "vjust", "rot", "sinRot", "cosRot",
                  "straightWidth", "straightHeight", "ascent", "descent",
                  "xShift", "yShift", "getFontsize", "refSize", "vpCex")
    if (!all(varNames %in% names(internals))) {
        stop(gettextf("'%s' must have the required variables",
                      'x[["internals"]]', domain = "R-sisal"), domain = NA)
    }
    oneGrob <- internals[["oneGrob"]]
    if (!resize && !identical(oneGrob, TRUE) && !identical(oneGrob, FALSE)) {
        stop("if 'resize' is FALSE, 'internals[[\"oneGrob\"]]' must be a logical flag")
    }
    rot <- internals[["rot"]]
    n <- length(rot)
    X <- internals[["x"]]
    Y <- internals[["y"]]
    hjust <- internals[["hjust"]]
    vjust <- internals[["vjust"]]
    sinRot <- internals[["sinRot"]]
    cosRot <- internals[["cosRot"]]
    straightWidth <- internals[["straightWidth"]]
    straightHeight <- internals[["straightHeight"]]
    ascent <- internals[["ascent"]]
    descent <- internals[["descent"]]
    sizingWidth <- internals[["sizingWidth"]]
    sizingHeight <- internals[["sizingHeight"]]
    widthOrig <- internals[["widthOrig"]]
    heightOrig <- internals[["heightOrig"]]
    refSize <- internals[["refSize"]]
    xShift <- internals[["xShift"]]
    yShift <- internals[["yShift"]]
    stopifnot(inherits(X, "unit"), inherits(Y, "unit"),
              is.numeric(hjust), is.numeric(vjust), is.numeric(sinRot),
              is.numeric(cosRot), is.numeric(straightWidth),
              is.numeric(straightHeight), is.numeric(ascent),
              is.numeric(descent), is.numeric(sizingWidth),
              is.numeric(sizingHeight), is.numeric(widthOrig),
              is.numeric(heightOrig), is.numeric(refSize), is.numeric(xShift),
              is.numeric(yShift), is.numeric(rot))
    stopifnot(length(X) == n, length(Y) == n,
              length(hjust) == n, length(vjust) == n, length(sinRot) == n,
              length(cosRot) == n, length(straightWidth) == n,
              length(straightHeight) == n, length(ascent) == n,
              length(descent) == n, length(sizingWidth) == 1,
              length(sizingHeight) == 1, length(widthOrig) == 1,
              length(heightOrig) == 1, length(refSize) == 1,
              length(xShift) == n, length(yShift) == n)
    stopifnot(is.finite(hjust), is.finite(vjust), is.finite(sinRot),
              is.finite(cosRot), is.finite(straightWidth),
              is.finite(straightHeight), is.finite(ascent),
              is.finite(descent), is.finite(sizingWidth),
              is.finite(sizingHeight), is.finite(widthOrig),
              is.finite(heightOrig), is.finite(refSize), is.finite(xShift),
              is.finite(yShift), is.finite(rot))

    charIdx <- internals[["charIdx"]]
    otherIdx <- internals[["otherIdx"]]
    if (!is.numeric(charIdx)) {
        stop(gettextf("'%s' must be numeric", 'internals[["charIdx"]]',
                      domain = "R-sisal"), domain = NA)
    }
    if (!is.numeric(otherIdx)) {
        stop(gettextf("'%s' must be numeric", 'internals[["otherIdx"]]',
                      domain = "R-sisal"), domain = NA)
    }
    if (length(charIdx) + length(otherIdx) != n) {
        stop(gettextf("'%s' and '%s' must total '%s' elements",
                      'internals[["charIdx"]]', 'internals[["otherIdx"]]',
                      'length(internals[["rot"]])', domain = "R-sisal"),
             domain = NA)
    }
    ## ## Testing that certain slots point to a function
    if (!is.function(internals[["getFontsize"]])) {
        stop(gettextf("'%s' must be a function", 'internals[["getFontsize"]]',
                      domain = "R-sisal"), domain = NA)
    }

    ## Testing getFontsize()
    fs <- internals[["getFontsize"]]()
    if (!is.numeric(fs) || length(fs) != 1 || !is.finite(fs)) {
        stop('internals[["getFontsize"]]() must return a fontsize')
    }

    ## Testing 'vpCex'
    vpCex <- internals[["vpCex"]]
    if (!resize &&
        !(is.numeric(vpCex) && length(vpCex) == 1 && is.finite(vpCex))) {
        stop("if 'resize' is FALSE, 'internals[[\"vpCex\"]]' must be a finite number")
    }

    ## Testing 'label'
    stopifnot(length(length(x[["label"]])) == 1)

    ## Testing 'x', 'y', 'xShift', 'yShift', 'width', 'height'
    X <- x[["x"]]
    Y <- x[["y"]]
    width <- x[["width"]]
    height <- x[["height"]]
    if (!inherits(X, "unit") || !inherits(Y, "unit") ||
        !inherits(width, "unit") || !inherits(height, "unit")) {
        stop("'x', 'y', 'width' and 'height' must be grid units")
    }
    if (length(width) != 1 || length(height) != 1) {
        stop("'width' and 'height' must have one element each")
    }

    ## Testing 'just', 'hjust', 'vjust'
    just <- x[["just"]]
    hjust <- x[["hjust"]]
    vjust <- x[["vjust"]]
    if (is.null(just)) {
        if (!is.numeric(hjust) || !is.finite(hjust)) {
            stop(gettextf("if '%s' is NULL, '%s' must be numeric and finite",
                          "just", "hjust", domain = "R-sisal"), domain = NA)
        }
        if (!is.numeric(vjust) || !is.finite(vjust)) {
            stop(gettextf("if '%s' is NULL, '%s' must be numeric and finite",
                          "just", "vjust", domain = "R-sisal"), domain = NA)
        }
    } else {
        if (!is.null(hjust) && (!is.numeric(hjust) || !is.finite(hjust))) {
            stop(gettextf("non-NULL '%s' must be numeric and finite", "hjust",
                          domain = "R-sisal"), domain = NA)
        }
        if (!is.null(vjust) && (!is.numeric(vjust) || !is.finite(vjust))) {
            stop(gettextf("non-NULL '%s' must be numeric and finite", "vjust",
                          domain = "R-sisal"), domain = NA)
        }
        if ((is.null(hjust) || is.null(vjust)) &&
            !(((is.numeric(just) && all(is.finite(just))) ||
               is.character(just)) && length(just) > 0 && length(just) < 3)) {
            stop("'just' must be a finite numeric or character vector of length 1 or 2")
        }
    }

    ## Testing 'sizingWidth', 'sizingHeight'
    sizingWidth <- x[["sizingWidth"]]
    sizingHeight <- x[["sizingHeight"]]
    if (!is.null(sizingWidth) &&
        (!is.numeric(sizingWidth) || !is.finite(sizingWidth) ||
         length(sizingWidth) != 1 || sizingWidth < 0)) {
        stop(gettextf("'%s' must be NULL or a non-negative number",
                      "sizingWidth", domain = "R-sisal"), domain = NA)
    }
    if (!is.null(sizingHeight) &&
        (!is.numeric(sizingHeight) || !is.finite(sizingHeight) ||
         length(sizingHeight) != 1 || sizingHeight < 0)) {
        stop(gettextf("'%s' must be NULL or a non-negative number",
                      "sizingHeight", domain = "R-sisal"), domain = NA)
    }

    ## Testing 'rot', 'rotJust', 'rotVjust', 'rotHjust'
    rot <- x[["rot"]]
    if (!is.numeric(rot) || length(rot) == 0 || !is.finite(rot)) {
        stop("'rot' must be a vector of finite numbers")
    }
    rotJust <- x[["rotJust"]]
    if (!is.logical(rotJust) || length(rotJust) == 0 ||
        !all(is.finite(rotJust))) {
        stop("'rotJust' must be a logical vector without any NA values")
    }
    rotVjust <- x[["rotVjust"]]
    if (!is.null(rotVjust) && !is.numeric(rotVjust)) {
        stop(gettextf("'%s' must be NULL or a numeric vector", "rotVjust",
                      domain = "R-sisal"), domain = NA)
    }
    rotHjust <- x[["rotHjust"]]
    if (!is.null(rotHjust) && !is.numeric(rotHjust)) {
        stop(gettextf("'%s' must be NULL or a numeric vector", "rotHjust",
                      domain = "R-sisal"), domain = NA)
    }

    ## Testing 'adjustJust'
    adjustJust <- x[["adjustJust"]]
    if (!(identical(adjustJust, TRUE) || identical(adjustJust, FALSE))) {
        stop(gettextf("'%s' must be a logical flag", "adjustJust",
                      domain = "R-sisal"), domain = NA)
    }

    x
}

editDetails.dynText <- function(x, specs, ...) {
    x2 <- NextMethod("editDetails", x)
    safeSlots <- c("x", "y", "name", "vp")
    recreate <- FALSE
    X <- x2[["x"]]
    Y <- x2[["y"]]
    n <- max(length(X), length(Y))
    internals <- x2[["internals"]]
    nOld <- length(internals[["x"]])
    slotNames <- names(specs)
    haveGp <- "gp" %in% slotNames
    haveResize <- "resize" %in% slotNames
    if ("internals" %in% slotNames) {
        warning("recreating grob because read-only slots were altered")
        recreate <- TRUE
    } else if (!all(slotNames %in% safeSlots) &&
               ((!x2[["resize"]] && x2[["adjustJust"]]) || n != nOld ||
                (x2[["resize"]] && (is.null(x2[["sizingWidth"]]) ||
                                    is.null(x2[["sizingHeight"]]))))) {
        recreate <- TRUE
    }
    if (recreate) {
        dynTextGrob(label = x2[["label"]], x = X, y = Y, width = x2[["width"]],
                    height = x2[["height"]], just = x2[["just"]],
                    hjust = x2[["hjust"]], vjust = x2[["vjust"]],
                    rot = x2[["rot"]], rotJust = x2[["rotJust"]],
                    rotVjust = x2[["rotVjust"]], rotHjust = x2[["rotHjust"]],
                    resize = x2[["resize"]], sizingWidth = x2[["sizingWidth"]],
                    sizingHeight = x2[["sizingHeight"]],
                    adjustJust = x2[["adjustJust"]],
                    name = x2[["name"]], gp = x2[["gp"]], vp = x2[["vp"]])
    } else {
        resize <- x2[["resize"]]
        X <- rep(x2[["x"]], length.out = n)
        Y <- rep(x2[["y"]], length.out = n)
        haveX <- "x" %in% slotNames
        haveY <- "y" %in% slotNames
        if (haveX) {
            internals[["x"]] <- X
        }
        if (haveY) {
            internals[["y"]] <- Y
        }
        if (!resize) {
            X <- X + unit(internals[["xShift"]], "inches")
            Y <- Y + unit(internals[["yShift"]], "inches")
        }
        editXY <- !resize && (haveX || haveY)
        if (any(c("just", "hjust", "vjust") %in% slotNames)) {
            internals[c("hjust", "vjust")] <-
                justHV(x2[["just"]], x2[["hjust"]], x2[["vjust"]], n)
            editJust <- TRUE
        } else {
            editJust <- FALSE
        }
        if ("rot" %in% slotNames) {
            rot <- x2[["rot"]]
            internals[["rot"]] <- rep_len(rot %% 360, n)
            rotRad <- rot * (pi / 180)
            internals[["sinRot"]] <- rep_len(sin(rotRad), n)
            internals[["cosRot"]] <- rep_len(cos(rotRad), n)
            editRot <- TRUE
        } else {
            editRot <- FALSE
        }
        rot2 <- internals[["rot"]]
        haveLabel <- "label" %in% slotNames
        gpList <- gparToList(x2[["gp"]])
        if (resize) {
            if (haveLabel) {
                label <- x2[["label"]]
                mfs <- plainLabelSize(label, n = n, rot = rot2,
                                      gpList = gpList, resize = resize,
                                      refSize = internals[["refSize"]])
                internals[["straightWidth"]] <- mfs[["width"]]
                internals[["straightHeight"]] <- mfs[["height"]]
                internals[["descent"]] <- mfs[["descent"]]
                internals[["ascent"]] <- mfs[["ascent"]]
                isCharacter <- mfs[["isCharacter"]]
                tmp <- setLabel.dynText(label, n, isCharacter, X, Y,
                                        internals[["hjust"]],
                                        internals[["vjust"]], rot2)
                x2 <- setChildren(x2, tmp[[1]])
                internals[["charIdx"]] <- tmp[[2]]
                internals[["otherIdx"]] <- tmp[[3]]
                internals[["oneGrob"]] <- tmp[[4]]
            } else if (internals[["oneGrob"]]) {
                if (editXY) {
                    x2 <- editGrob(x2, "one", x = X, y = Y)
                }
                if (editJust) {
                    hjust <- internals[["hjust"]]
                    vjust <- internals[["vjust"]]
                    x2 <- editGrob(x2, "one", hjust = hjust, vjust = vjust)
                }
                if (editRot) {
                    x2 <- editGrob(x2, "one", rot = rot2)
                }
            } else if (editXY || editJust || editRot) {
                charIdx <- internals[["charIdx"]]
                otherIdx <- internals[["otherIdx"]]
                if (editXY) {
                    x2 <- editGrob(x2, "char", x = X[charIdx], y = Y[charIdx])
                    x2 <- editGrob(x2, "other",
                                   x = X[otherIdx], y = Y[otherIdx])
                }
                if (editJust) {
                    hjust <- internals[["hjust"]]
                    vjust <- internals[["vjust"]]
                    x2 <- editGrob(x2, "char", hjust = hjust[charIdx],
                                   vjust = vjust[charIdx])
                    x2 <- editGrob(x2, "other", hjust = hjust[otherIdx],
                                   vjust = vjust[otherIdx])
                }
                if (editRot) {
                    x2 <- editGrob(x2, "char", rot = rot2[charIdx])
                    x2 <- editGrob(x2, "other", rot = rot2[otherIdx])
                }
            }
        }
        if ("sizingWidth" %in% slotNames) {
            newValue <- x2[["sizingWidth"]]
            if (is.null(newValue)) {
                internals[["sizingWidth"]] <- internals[["widthOrig"]]
            } else {
                internals[["sizingWidth"]] <- newValue
            }
        }
        if ("sizingHeight" %in% slotNames) {
            newValue <- x2[["sizingHeight"]]
            if (is.null(newValue)) {
                internals[["sizingHeight"]] <- internals[["heightOrig"]]
            } else {
                internals[["sizingHeight"]] <- newValue
            }
        }
        if (resize) {
            gpList[["fontsize"]] <- NULL
            gpList[["cex"]] <- 1
        }
        if (haveResize || haveLabel || haveGp ||
            any(c("width", "height", "sizingWidth",
                  "sizingHeight") %in% slotNames)) {
            internals[["getFontsize"]] <-
                makeGetFontsize.dynText(x2[["label"]], gpList,
                                        internals[["sizingWidth"]],
                                        internals[["sizingHeight"]],
                                        internals[["refSize"]], n,
                                        x2[["width"]], x2[["height"]], resize)
        }
        if (!resize) {
            if (haveResize || haveGp) {
                currentGP <- gparToList(get.gpar())
                vpCex <- currentGP[["cex"]]
                if (is.null(vpCex)) {
                    vpCex <- 1
                } else {
                    currentGP[["cex"]] <- 1
                }
                internals[["vpCex"]] <- vpCex
            }
            if (haveGp) {
                currentGP[names(gpList)] <- gpList
                gpList <- gparToList(currentGP)
                x2[["gp"]] <- do.call("gpar", gpList)
            }
        } else if (haveGp) {
            x2[["gp"]] <- do.call("gpar", gpList)
        }
        x2[["internals"]] <- internals
        x2
    }
}
