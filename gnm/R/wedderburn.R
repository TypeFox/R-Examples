#  Modification of binomial from the stats package for R.
#
#  Copyright (C) 1995-2005 The R Core Team
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

"wedderburn" <-
    function (link = "logit")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link")
            linktemp <- eval(link)
    }
    if (any(linktemp == c("logit", "probit", "cloglog")))
        stats <- make.link(linktemp)
    else stop(paste(linktemp,
                    "link not available for wedderburn quasi-family;",
                    "available links are",
                    "\"logit\", \"probit\" and \"cloglog\""))
    variance <- function(mu)  mu^2 * (1-mu)^2
    validmu <- function(mu) {
        all(mu > 0) && all(mu < 1)}
    dev.resids <- function(y, mu, wt){
        eps <-  0.0005
        2 * wt * (y/mu + (1 - y)/(1 - mu) - 2 +
                  (2 * y - 1) * log((y + eps)*(1 - mu)/((1- y + eps) * mu)))
    }
    aic <- function(y, n, mu, wt, dev) NA
    initialize <- expression({
        if (any(y < 0 | y > 1)) stop(paste(
                   "Values for the wedderburn family must be in [0,1]"))
        n <- rep.int(1, nobs)
        mustart <- (y + 0.1)/1.2
    })
    structure(list(family = "wedderburn",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta),
              class = "family")
}
