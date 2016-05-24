#  Copyright (C) 2005, 2008, 2013 Heather Turner
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

residuals.gnm <- function(object, type = "deviance", ...) {
    if (type == "partial")
        stop("type = \"partial\" not implemented for gnm objects.")
    else res <- NextMethod("residuals")

    if (!is.null(object$table.attr))
        attributes(res) <- object$table.attr

    res
}
