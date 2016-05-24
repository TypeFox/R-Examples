#  Copyright (C) 2005, 2006, 2008 Heather Turner
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

Mult <- function(..., inst = NULL){
    if ("multiplicity" %in% names(match.call()[-1]))
        stop("multiplicity argument of Mult has been replaced by",
             "\"inst\" argument.")
    dots <- match.call(expand.dots = FALSE)[["..."]]
    list(predictors = dots,
         term = function(predLabels, ...) {
             paste("(", paste(predLabels, collapse = ")*("), ")", sep = "")
         },
         call = as.expression(match.call()),
         match = seq(dots))
}
class(Mult) <- "nonlin"
