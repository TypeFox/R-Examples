#  Modification of model.frame.glm from the stats package for R.
#
#  Copyright (C) 1995-2005 The R Core Team
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

model.frame.gnm <- function (formula, ...)
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots),
        0)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1]] <- as.name("gnm")
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env))
            env <- parent.frame()
        eval(fcall, env)
    }
    else formula$model
}



