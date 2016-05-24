#  Modification of confint.glm from the MASS package for R.
#
#  Copyright (C) 1994-2006 W. N. Venables and B. D. Ripley
#  Copyright (C) 2006 Heather Turner
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

confint.gnm <- function (object, parm = ofInterest(object), level = 0.95,
                         trace = FALSE, ...)
{
    pnames <- names(coef(object))
    if (is.null(parm))
        parm <- seq(along = pnames)
    else if (is.character(parm))
        parm <- match(parm, pnames, nomatch = 0)
    cat("Waiting for profiling to be done...\n")
    flush.console()
    object <- profile(object, which = parm, alpha = 1 - level, trace = trace)
    confint(object, level = level, ...)
}
