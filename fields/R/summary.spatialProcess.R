# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"summary.spatialProcess" <- function(object, digits = 4, ...) {
    x <- object
 
    summary <- list(call = x$call, num.observation = length(x$residuals), 
        enp = x$eff.df, nt = x$nt, df.drift = sum(x$ind.drift), 
        res.quantile = quantile(x$residuals, seq(0, 1, 0.25)), 
        sigma.MLE = x$shat.MLE, theta.MLE = x$args$theta, 
        shat.GCV = x$shat.GCV, rho.MLE = x$rho.MLE, 
        m = x$m, lambda = x$lambda, cost = x$cost, rho = x$rho, 
        sigma2 = x$sigma2, num.uniq = length(x$yM), knot.model = x$knot.model, 
        np = x$np, method = x$method,  shat.pure.error = x$shat.pure.error, 
        args = x$args)
    class(summary) <- "summary.spatialProcess"
    summary$covariance <- cor(x$fitted.values * sqrt(x$weights), 
        (x$y) * sqrt(x$weights))^2
    hold <- (sum((x$y - mean(x$y))^2) - sum(x$residuals^2))/(sum((x$y - 
        mean(x$y))^2))
    summary$adjr2 <- 1 - ((length(x$residuals) - 1)/(length(x$residuals) - 
        x$eff.df)) * (1 - hold)
    summary$digits <- digits
    summary$cov.function <- as.character(x$cov.function.name)
    summary$correlation.model <- x$correlation.model
    summary$sum.gcv.lambda <- x$lambda.est[6,]
    summary
}
