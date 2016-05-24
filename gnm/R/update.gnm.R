#  Modification of update.default from the stats package for R.
#
#  Copyright (C) 1995-2010 The R Core Team
#  Copyright (C) 2010, 2012 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

update.gnm <- function (object, formula., ..., evaluate = TRUE)
{
    call <- object$call
    if (is.null(call))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) {
        ## update.formula reorders nonlin terms as lin (main effects)
        ## therefore use substitute to keep order
        formula. <- as.formula(formula.)
        rhs <- formula.[[length(formula.)]]
        rhs <- do.call(substitute,
                       list(rhs, env = list("." = object$formula[[3]])))
        if (length(formula.) == 3) {
            lhs <- formula.[[2]]
            lhs <- do.call(substitute,
                           list(lhs, env = list("." = object$formula[[2]])))
            call$formula <- call("~", lhs, rhs)
        } else call$formula <- call("~", object$formula[[2]], rhs)
        call$formula <- formula(terms.formula(call$formula, simplify = TRUE,
                                              keep.order = TRUE))
    }
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}

