### File R/plotSelected.R
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

## Given a character vector x, extracts the integer part of matching
## member strings.  A string matches if only one integer part is
## present.  Other strings are left unchanged.
intPart <- function(x, tform = function(x) x, ...) {
    haveNumPart <- grep("^\\D*[\\+-]?\\d+\\D*$", x)
    y <- x
    numPart <- tform(as.integer(sub("^\\D*([\\+-]?\\d+)\\D*$", "\\1",
                                    x[haveNumPart])))
    y[haveNumPart] <- formatC(x = numPart, ...)
    y
}

.plotSelected <- function(x, useAllNames = TRUE,
                          pickIntPart = FALSE,
                          intTransform = function(x) x,
                          formatCArgs = list(),
                          xLabels = 1, yLabels = NULL,
                          L.f.color = "black", L.v.color = "grey50",
                          other.color = "white",
                          naFill = other.color, naStripes = L.v.color,
                          selectedLabels = TRUE, otherLabels = FALSE,
                          labelPar = gpar(fontface = 1, fontsize = 20,
                          cex = 0.35),
                          nestedPar = gpar(fontface = 3),
                          ranking = c("pairwise", "nested"),
                          tableArgs = list(), ...) {
    ## If nested ranking exists AND
    ## "nested" %in% ranking:
    ## * use labelPar + nestedPar (nestedPar overwrites labelPar)
    ## * if both ranking methods are requested in argument 'ranking':
    ##   - let A be the ranking referred to in the first element of
    ##     argument 'ranking' and B be the other ranking
    ##   - if A[k] != B[k], use label sprintf("%.0f (%.0f)", A[k], B[k])
    ##   - if A[k] == B[k], use label sprintf("%.0f", A[k])
    ## Otherwise:
    ## * use labelPar only
    ## * use pairwise ranking only
    stopifnot(identical(useAllNames, TRUE) || identical(useAllNames, FALSE),
              identical(selectedLabels, TRUE) ||
              identical(selectedLabels, FALSE),
              identical(otherLabels, TRUE) || identical(otherLabels, FALSE),
              is.logical(pickIntPart), length(pickIntPart) > 0,
              pickIntPart %in% c(TRUE, FALSE),
              is.function(intTransform), is.list(formatCArgs),
              length(formatCArgs) == 0 || !is.null(names(formatCArgs)),
              inherits(labelPar, "gpar"), inherits(nestedPar, "gpar"),
              is.list(tableArgs))
    ranking2 <- unique(match.arg(ranking, several.ok = TRUE))
    if (inherits(x, "sisal")) {
        x2 <- list(x)
    } else if (is.list(x)) {
        foo <- vapply(x, inherits, TRUE, what = "sisal")
        if (any(!foo)) {
            stop("'x' must only contain \"sisal\" objects")
        }
        x2 <- x
    } else {
        stop("'x' must be a \"sisal\" object or a list")
    }
    nRows <- length(x2)
    if (nRows == 0) {
        stop("no \"sisal\" objects given")
    }
    nPick <- length(pickIntPart)
    if (useAllNames) {
        xLabels2 <- character(0)
        inputMapping <- vector(mode = "list", length = nRows)
        for (k in seq_len(nRows)) {
            theseNames <- x2[[k]][["var.names"]]
            if (is.null(theseNames)) {
                stop(gettextf("no input names in sisal object number %.0f",
                              k, domain="R-sisal"), domain=NA)
            }
            if (pickIntPart[(k - 1) %% nPick + 1]) {
                fArgs <- formatCArgs
                fArgs[["x"]] <- theseNames
                fArgs[["tform"]] <- intTransform
                theseNames <- do.call("intPart", fArgs)
            }
            xLabels2 <- c(xLabels2,
                          theseNames[!theseNames %in% xLabels2])
            inputMapping[[k]] <- match(theseNames, xLabels2)
        }
        nCols <- length(xLabels2)
        if ((is.character(xLabels) || is.list(xLabels))) {
            if (all(xLabels2 %in% names(xLabels))) {
                matches <- match(xLabels2, names(xLabels))
                newOrder <- sort.list(matches)
                xLabels2 <- xLabels[matches[newOrder]]
                for (k in seq_len(nRows)) {
                    inputMapping[[k]] <- match(inputMapping[[k]], newOrder)
                }
            } else {
                warning("'xLabels' does not have names")
                xLabels2 <- xLabels
                if (nCols > length(xLabels2)) {
                    warning("x axis labels not available for all columns")
                    xLabels2 <- c(xLabels2,
                                  rep.int(NA_character_,
                                          nCols - length(xLabels2)))
                } else if (nCols < length(xLabels2)) {
                    warning(gettextf("extra elements of '%s' ignored",
                                     "xLabels", domain = "R-sisal"),
                            domain = NA)
                    xLabels2 <- xLabels2[seq_len(nCols)]
                }
            }
        }
    } else if (is.numeric(xLabels)) {
        if (!(length(xLabels) == 1 && is.finite(xLabels) &&
              round(xLabels) == xLabels && xLabels >= 1 && xLabels <= nRows)) {
            stop("numeric 'xLabels' must be an index to 'x'")
        }
        xLabels2 <- x2[[xLabels]][["var.names"]]
        if (is.null(xLabels2)) {
            stop(gettextf("no input names in sisal object number %.0f",
                          xLabels, domain="R-sisal"), domain=NA)
        }
        if (pickIntPart[(xLabels - 1) %% nPick + 1]) {
            fArgs <- formatCArgs
            fArgs[["x"]] <- xLabels2
            fArgs[["tform"]] <- intTransform
            xLabels2 <- do.call("intPart", fArgs)
        }
        nInputs <- vapply(x2, function (x) x[["d"]], 1)
        inputMapping <- lapply(nInputs, seq_len)
        nCols <- max(nInputs)
        if (any(nInputs != nCols)) {
            warning("sisal objects have different numbers of inputs")
            if (nCols > length(xLabels2)) {
                warning("x axis labels not available for all columns")
                xLabels2 <- c(xLabels2,
                              rep.int(NA_character_, nCols - length(xLabels2)))
            }
        }
    } else if (is.character(xLabels) || is.list(xLabels)) {
        xLabels2 <- xLabels
        nInputs <- vapply(x2, function (x) x[["d"]], 1)
        inputMapping <- lapply(nInputs, seq_len)
        nCols <- max(nInputs)
        if (any(nInputs != nCols)) {
            warning("sisal objects have different numbers of inputs")
        }
        if (nCols > length(xLabels2)) {
            warning("x axis labels not available for all columns")
            xLabels2 <- c(xLabels2,
                          rep.int(NA_character_, nCols - length(xLabels2)))
        } else if (nCols < length(xLabels2)) {
            warning(gettextf("extra elements of '%s' ignored",
                             "xLabels", domain = "R-sisal"),
                    domain = NA)
            xLabels2 <- xLabels2[seq_len(nCols)]
        }
    } else {
        stop("if 'not(useAllNames)', 'xLabels' must be numeric, character or list")
    }
    if (nCols == 0) {
        stop("no inputs in sisal objects (this should not happen)")
    }
    if (is.null(yLabels)) {
        yLabels2 <- NULL
    } else if (is.character(yLabels) || is.list(yLabels)) {
        yLabels2 <- yLabels
        if (nRows > length(yLabels2)) {
            warning("'yLabels' not available for all sisal objects")
            yLabels2 <- c(yLabels2, rep.int("", nRows - length(yLabels2)))
        } else if (nRows < length(yLabels2)) {
            warning(gettextf("extra elements of '%s' ignored",
                             "yLabels", domain = "R-sisal"),
                    domain = NA)
            yLabels2 <- yLabels2[seq_len(nRows)]
        }
    } else {
        stop("'yLabels' must be NULL, character or list")
    }
    ## Now we have nRows, nCols, xLabels2, yLabels2 and inputMapping
    nCells <- nRows * nCols
    bg <- rep.int(NA_character_, nCells)
    stripeCol <- rep.int(NA_character_, nCells)
    labels <- rep.int(NA_character_, nCells)
    labelPars <- vector(mode="list", length=nCells)
    labelPars[] <- list(gpar())
    nestedPar2 <- gparToList(labelPar)
    bPar <- gparToList(nestedPar)
    for (name in names(bPar)) {
        nestedPar2[[name]] <- bPar[[name]]
    }
    nestedPar2 <- do.call("gpar", nestedPar2)
    nestedIdx <- which(ranking2 == "nested")
    pairIdx <- which(ranking2 == "pairwise")
    haveNested <- length(nestedIdx) > 0
    havePair <- length(pairIdx) > 0
    if (haveNested && havePair) {
        nestedFirst <- nestedIdx < pairIdx
    }
    for (k in seq_len(nRows)) {
        thisSisal <- x2[[k]]
        nInputs <- thisSisal[["d"]]
        li <- (inputMapping[[k]] - 1) * nRows + k
        L.v <- thisSisal[["L.v"]]
        L.f <- thisSisal[["L.f"]]
        if (haveNested && length(thisSisal[["path.length"]]) == 1) {
            if (havePair) {
                if (nestedFirst) {
                    rankingA <- thisSisal[["nested.rank"]]
                    rankingB <- thisSisal[["pairwise.rank"]]
                } else {
                    rankingA <- thisSisal[["pairwise.rank"]]
                    rankingB <- thisSisal[["nested.rank"]]
                }
                diff.flag <- rankingA != rankingB
                rankChar <- character(nInputs)
                rankChar[diff.flag] <- sprintf("%.0f (%.0f)",
                                               rankingA[diff.flag],
                                               rankingB[diff.flag])
                rankChar[!diff.flag] <- sprintf("%.0f", rankingA[!diff.flag])
            } else {
                rankChar <- sprintf("%.0f", thisSisal[["nested.rank"]])
            }
            labelPars[li] <- list(nestedPar2)
        } else {
            rankChar <- sprintf("%.0f", thisSisal[["pairwise.rank"]])
            labelPars[li] <- list(labelPar)
        }
        bg[li] <- other.color
        L.vf <- intersect(L.v, L.f)
        L.vOnly <- setdiff(L.v, L.vf)
        L.fOnly <- setdiff(L.f, L.vf)
        L.vOrF <- union(L.f, L.v)
        bg[li[L.vOnly]] <- L.v.color
        bg[li[L.fOnly]] <- L.f.color
        bg[li[L.vf]] <- L.f.color
        stripeCol[li[L.vf]] <- L.v.color
        if (selectedLabels) {
            labels[li[L.vOrF]] <- rankChar[L.vOrF]
        }
        if (otherLabels) {
            L.others <- rep.int(TRUE, nInputs)
            L.others[L.vOrF] <- FALSE
            labels[li[L.others]] <- rankChar[L.others]
        }
    }
    cl <- as.call(c(as.name("sisalTable"), c(list(...), tableArgs)))
    fArgs <- as.list(match.call(sisalTable, cl))[-1L]
    fArgs[["nRows"]] <- nRows
    fArgs[["nCols"]] <- nCols
    fArgs[["bg"]] <- bg
    fArgs[["stripeCol"]] <- stripeCol
    fArgs[["labels"]] <- labels
    fArgs[["labelPars"]] <- labelPars
    fArgs[["xAxisLabels"]] <- xLabels2
    fArgs[["yAxisLabels"]] <- yLabels2
    fArgs[["naFill"]] <- naFill
    fArgs[["naStripes"]] <- naStripes
    do.call("sisalTable", fArgs)
}

setMethodS3("plotSelected", "sisal", .plotSelected)
setMethodS3("plotSelected", "list", function(x, ...) {
    ## Normal function call, no method dispatch.  Selects the function
    ## based on the class of the first element of the list.  Thus,
    ## possible additional "plotSelected" methods should work:
    ## * with one object or a list of objects belonging to the same class
    ## * through method dispatch or when called directly (i.e. no
    ##   reliance on the special environment set up by the method
    ##   dispatch system)
    getDispatchMethodS3("plotSelected", class(x[[1]]))(x, ...)
})
