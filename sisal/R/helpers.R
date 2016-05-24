### File R/helpers.R
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

justHV <- function(just, hjust = NULL, vjust = NULL, n = 1) {
    if (is.null(hjust)) {
        if (is.numeric(just)) {
            hjust2 <- just[1]
        } else if (just[1] == "left") {
            hjust2 <- 0
        } else if (just[1] == "right") {
            hjust2 <- 1
        } else {
            hjust2 <- 0.5
        }
        Hjust <- rep.int(hjust2, n)
    } else {
        Hjust <- rep_len(hjust, n)
    }
    if (is.null(vjust)) {
        vjust2 <- 0.5
        if (is.numeric(just)) {
            if (length(just) > 1) {
                vjust2 <- just[2]
            }
        } else if (length(just) == 1) {
            if (just == "top") {
                vjust2 <- 1
            } else if (just == "bottom") {
                vjust2 <- 0
            }
        } else if (just[2] == "top") {
            vjust2 <- 1
        } else if (just[2] == "bottom") {
            vjust2 <- 0
        }
        Vjust <- rep.int(vjust2, n)
    } else {
        Vjust <- rep_len(vjust, n)
    }
    list(Hjust, Vjust)
}

vars.to.name <- function(vars.boolean, empty = FALSE) {
    if (empty) {
        "0"
    } else {
        paste0(as.character(which(vars.boolean)), collapse=".")
    }
}

name.to.vars <- function(name, d) {
    res <- rep.int(FALSE, d)
    if (name != "0") {
        res[as.numeric(strsplit(name, ".", fixed=TRUE)[[1]])] <- TRUE
    }
    res
}

name.to.vars.idx <- function(name) {
    if (name != "0") {
        as.numeric(strsplit(name, ".", fixed=TRUE)[[1]])
    } else {
        numeric(0)
    }
}

name.to.vars.char <- function(name) {
    if (name != "0") {
        strsplit(name, ".", fixed=TRUE)[[1]]
    } else {
        character(0)
    }
}

splitcat <- function(...) {
    for(splitString in strsplit(unlist(list(...)), " ")) {
        cat(splitString, fill=TRUE)
    }
}

digestVector <- function(x) {
    toChar <- function(x) {
        if (is.double(x)) {
            xChar <- formatC(x, digits=17, width=-1, format="g")
        } else if (is.integer(x)) {
            xChar <- paste0(formatC(x, format="d"), "L")
        } else if (is.list(x)) {
            xChar <- c("(", vapply(x, toChar, character(1)), ")")
        } else if (is.character(x)) {
            xChar <- paste0('"', x, '"')
        } else {
            xChar <- as.character(x)
        }
        paste0(xChar, collapse=" ")
    }
    digest(toChar(x), algo="md5", serialize=FALSE)
}

## Returns TRUE for atomic vectors and calls to a limited set of
## functions: sprintf, gettextf, gettext and ngettext in the base
## namespace.
safeToEval <- function(object, depth=0) {
    if (depth >= 10) {
        FALSE # reasonable limit for recursion depth
    } else if (is.atomic(object)) {
        TRUE
    } else if (is.call(object)) {
        GOOD.FUNCS <- c("sprintf", "gettextf", "gettext", "ngettext")
        foo <- as.list(object)
        foo.1 <- foo[[1]]
        n.foo <- length(foo)
        if (n.foo < 2) {
            FALSE
        } else if (is.name(foo.1)) {
            if (as.character(foo.1) %in% GOOD.FUNCS) {
                all(vapply(foo[2:n.foo], safeToEval, TRUE, depth=depth + 1))
            } else {
                FALSE
            }
        } else if (is.call(foo.1)) {
            bar <- as.list(foo.1)
            n.bar <- length(bar)
            if (n.bar != 3) {
                FALSE
            } else if (is.name(bar[[1]]) &&
                       as.character(bar[[1]]) %in% c("::", ":::") &&
                       (is.name(bar[[2]]) || is.character(bar[[2]])) &&
                       identical(as.character(bar[[2]]), "base") &&
                       (is.name(bar[[3]]) ||
                        (is.character(bar[[3]]) && length(bar[[3]])==1)) &&
                       as.character(bar[[3]]) %in% GOOD.FUNCS) {
                all(vapply(foo[2:n.foo], safeToEval, TRUE, depth=depth + 1))
            } else {
                FALSE
            }
        } else {
            FALSE
        }
    } else {
        FALSE
    }
}

## There seems to be no documented way of modifying a gpar object. The
## user should convert the gpar to a list, then modify it and create a
## new gpar object. This function removes "font" if both "font" and
## "fontface" are specified, Unnamed members are also removed.
gparToList <- function(gp) {
    gpList <- unclass(gp)
    theNames <- names(gpList)
    if (is.null(theNames)) {
        list()
    } else {
        if (all(c("font", "fontface") %in% theNames)) {
            gpList["font"] <- NULL
            theNames <- names(gpList)
        }
        gpList[!is.na(theNames) & nchar(theNames) > 0]
    }
}

## This is like gparToList(), but also edits the parameters (name =
## value) given in ...
editGpar <- function(gp, ...) {
    gpList <- gparToList(gp)
    args <- list(...)
    argNames <- names(args)
    nArgs <- length(args)
    argHasName <- logical(nArgs)
    if (!is.null(argNames)) {
        argHasName[!is.na(argNames) & nchar(argNames) > 0] <- TRUE
        gpList[argNames[argHasName]] <- args[argHasName]
        do.call("gpar", gpList)
    } else {
        gp
    }
}

## Determines if range of vector is zero (with tolerance) or NA.
## Based on a function by Hadley Wickham,
## http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
zeroRange <- function(x, mean.x = mean(x, na.rm=TRUE),
                      range.x = suppressWarnings(range(x, na.rm=TRUE)),
                      tol = .Machine[["double.eps"]] ^ 0.5) {
    xRange <- range.x / mean.x
    isTRUE(all.equal(xRange[1], xRange[2], tolerance = tol))
}
