#  File permute/R/update.how.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
# Modifications by Gavin L. Simpson
#
# Copyright (C) 2013 Gavin L. Simpson
#
# Modifcations made:
#  1) Remove `formula.` argument and processing thereof.
#  2) Evaluation is forced to the global environment.
#  3) (minor) Took the definition of `call` out of the `if ()` statement
#     for clarity/style issues.
#  4) Added this modification section to the copyright/licence header.
#  5) Added code to preserve some components of the original object.

`update.how` <- function (object, ..., evaluate = TRUE) {
    call <- getCall(object)
    if (is.null(call))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...

    ## preserve or update block and/or plot  names
    bname <- if ("blocks" %in% names(extras)) {
        deparse(substitute(extras[["blocks"]]))
    } else {
        object$blocks.name
    }
    pname <- if ("plots" %in% names(extras)) {
        dots <- list(...)
        dots$plots$plots.name
    } else {
        object$plots$plots.name
    }

    ## process remaining ... args
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        ## do these individually to allow NULL to remove entries.
        for (a in names(extras)[existing])
            call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }

    ## probably want to evaluate hence default is TRUE
    if (evaluate) {
        out <- eval(call, parent.frame())
        ## add back in the chars we discovered earlier
        out$blocks.name <- bname
        out$plots$plots.name <- pname
    } else {
        out <- call
    }

    out
}
