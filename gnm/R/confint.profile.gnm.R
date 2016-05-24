#  Modification of confint.profile.glm from the MASS package for R.
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

confint.profile.gnm <- function (object, parm = names(object),
                                 level = 0.95, ...)  {
    of <- attr(object, "original.fit")
    pnames <- names(coef(of))
    if (is.numeric(parm))
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(round(100 * a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
    cutoff <- qnorm(a)
    std.err <- attr(object, "summary")$coefficients[, "Std. Error"]
    parm <- parm[!is.na(std.err)[parm]]
    for (pm in parm) {
        pro <- object[[pm]]
        if (is.matrix(pro[, "par.vals"]))
            sp <- spline(x = pro[, "par.vals"][, pm], y = pro[,
                1])
        else sp <- spline(x = pro[, "par.vals"], y = pro[, 1])
        print(pro[, "par.vals"][, pm])
        print(pro[,1])
        est <- approx(sp$y, sp$x, xout = cutoff)$y
        ci[pm, ] <- ifelse(is.na(est) & attr(pro, "asymptote"),
                           c(-Inf, Inf), est)
    }
    drop(ci)
}
