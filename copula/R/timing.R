## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### Timing for the implemented nested Archimedean copulas

##' Computes user times for the admissible parameter combinations provided by "taus"
##'
##' @title Timing frailties
##' @param n number of variates to be generated
##' @param family the (nested) Archimedean family to be timed
##' @param taus the sequence of Kendall's tau to be tested
##' @param digits number of digits for the output
##' @param verbose print current state of the timing during timing
##' @return a (tau_0 x tau_1)-matrix with first column indicating the user run
##'         times for V0 and the other cells the run time for V01 corresponding
##'         to two given taus among "taus" based on the generated V0's
##' @author Marius Hofert
nacFrail.time <- function(n, family, taus, digits=3, verbose=FALSE)
{
    ## setup
    f <- function(x) formatC(x, digits=digits, width = 8)
    mTime <- function(x) 1000 * system.time(x)[1] # measuring milliseconds
    l <- length(taus)
    f.taus <- format(taus, digits=digits)
    res <- matrix(,nrow=l,ncol=l,
		  ## use taus as row and column headers:
		  dimnames =
		  list('outer tau' = f.taus,
		       '   inner tau' = c(" ", f.taus[-1])))
    cop <- getAcop(family)
    thetas <- cop@iTau(taus)

    ## timing (based on user time)
    for(i in seq_along(thetas)) { # run over all theta0
        ## run times for V0 go into the first column:
        res[i,1] <- mTime(V0 <- cop@V0(n,thetas[i]))
        if(verbose) cat("V0:  tau_0 = ",f.taus[i],
                        "; time = ", f(res[i,1]), " ms\n",sep="")
        if(i < l) for(j in (i+1):l) {   # run over all theta1
            res[i,j] <- mTime(V01 <- cop@V01(V0,thetas[i], thetas[j]))
            if(verbose) cat("  V01: tau_0 = ",f.taus[i],", tau_1 = ", f.taus[j],
                            "; time = ",f(res[i,j])," ms\n", sep="")
        }
    }
    res
}
