#  Modification of cooks.distance.glm from the stats package for R.
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

cooks.distance.gnm <- function(model, hat = hatvalues(model), dispersion =
                               attr(vcov(model), "dispersion"), ...){
    p <- model$rank
    res <- na.omit(residuals(model, type =
                             "pearson"))[model$prior.weights != 0]
    res <- naresid(model$na.action, res)
    res <- (res/(1 - hat))^2 * hat/(dispersion * p)
    res[is.infinite(res)] <- NaN
    res
}
