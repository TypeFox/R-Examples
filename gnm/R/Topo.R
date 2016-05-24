#  Copyright (C) 2005 David Firth
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

"Topo" <-
    function (..., spec = NULL)
{
    if (is.null(spec)) stop("No spec given")
    dots <- list(...)
    factorLengths <- sapply(dots, length)
    lengthsEqual <- {if (length(factorLengths) == 1) TRUE else
                     sd(factorLengths) == 0}
    if (!lengthsEqual) stop("Factors have different lengths")
    specDim <- if (is.vector(spec)) length(spec) else dim(spec)
    dots <- lapply(dots, as.factor)
    facMat <- cbind(...)
    spec.ok <- identical(sapply(dots, nlevels), specDim)
    if (!spec.ok) stop(
            "Dimensions of spec do not match the factor arguments")
    return(as.factor(spec[facMat]))
}
