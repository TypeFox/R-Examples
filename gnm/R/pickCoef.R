#  Copyright (C) 2006, 2010, 2012, 2013 Heather Turner
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

pickCoef <- function(object, pattern = NULL, value = FALSE, ...){
    coefs <- names(coef(object))
    if (is.null(coefs))
        stop("Coefficient names cannot be extracted from 'object'")
    if (missing(pattern)) {
        default <- list(setlabels = "Selected coefficients",
                        title = "Select coefficients of interest",
                        items.label = "Model coefficients:",
                        return.indices = TRUE, edit.setlabels = FALSE,
                        warningText =  "No subset of coefficients selected")
        dots <- list(...)
        dotArgs <- match(names(default), names(dots))
        allArgs <- c(list(coefs), dots, default[is.na(dotArgs)])
        selection <- unname(unlist(do.call(pickFrom, allArgs)))
    }
    else {
        selection <- grep(pattern, coefs, value = FALSE, ...)
    }

    if (!length(selection))
        selection <- NULL
    else if (!value)
        names(selection) <- coefs[selection]
    else
        selection <- parameters(object)[selection]
    selection
}
