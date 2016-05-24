#  Copyright (C) 2005, 2006 Heather Turner
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

model.matrix.gnm <- function(object, coef = NULL, ...) {
    if (!"x" %in% names(object) || !is.null(coef)) {
        xcall <- object$call
        xcall$method <- "model.matrix"
        xcall$constrain <- object$constrain
        xcall$constrainTo <- object$constrainTo
        xcall$data <- model.frame(object)
        xcall[c("weights", "offset")] <- NULL
        xcall$verbose <- FALSE
        if (!is.null(coef))
            xcall$start <- coef
        else
            xcall$start <- coef(object)
        eval(xcall)
    }
    else object[[match("x", names(object))]]
}
