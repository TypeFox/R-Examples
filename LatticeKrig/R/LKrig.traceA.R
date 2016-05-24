# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

LKrig.traceA <- function(Mc, wX, wU, lambda, weights, NtrA, iseed = NA) {
# do not disrupt an external use of the random number generator    	
    if (exists(".Random.seed", 1)) {
        save.seed <- .Random.seed
    }
    else {
        save.seed <- NA
    }
    n <- length(weights)
    # set  new seed to use for Monte Carlo estimate of trace A(lambda)
    if (!is.na(iseed)) {
        set.seed(iseed)
    }
    # generate N(0,1)
    wEy <- matrix(rnorm(NtrA * n), n, NtrA) * sqrt(weights)
    # restore  original random seed
    if (!is.na(iseed) & !is.na(save.seed[1])) {
        assign(".Random.seed", save.seed, pos = 1)
    }
    #
    out3 <- LKrig.coef(Mc, wX, wU, wEy, lambda) 
    wEyhat <- (wX %*% out3$c.coef)
    if( !is.null(wU) ){
        wEyhat <- wEyhat + wU %*% out3$d.coef
    }
    trA.info <- colSums((wEy * wEyhat)/weights)
    trA.est <- mean(trA.info)
    trA.SE <- sqrt(var(trA.info)/length(trA.info))
    list(trA.est = trA.est, eff.df = trA.est,  trA.SE = trA.SE)
}
