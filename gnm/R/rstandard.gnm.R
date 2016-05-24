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

rstandard.gnm <- function(model, ...) {
    so <- summary(model)
    res <- na.omit(so$deviance.resid[model$prior.weights != 0])
    res <- naresid(model$na.action, res)
    res <- res/sqrt(so$dispersion * (1 - hatvalues(model)))
    res[is.infinite(res)] <- NaN
    if (!is.null(model$table.attr))
        attributes(res) <- model$table.attr
    res
}
