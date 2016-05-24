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


###  Generating Exponentially Tilted Stable Random Vars (r ET stable)
###  ================================================================
###  Experiments with retstable*() versions
###
###  More (computer intensive) experiments are in
###  ../demo/retstable.R
###  ~~~~~~~~~~~~~~~~~~~

require(copula)
source(system.file("Rsource", "utils.R", package="copula"))
##--> tryCatch.W.E(), canGet()

## it works for 0-length  V0 as well:
.N <- numeric(0) ; stopifnot(identical(.N, retstable(1/4, .N)))

## This is from "next version of Matrix" test-tools-1.R:
showSys.time <- function(expr) {
    ## prepend 'Time' for R CMD Rdiff
    st <- system.time(expr)
    writeLines(paste("Time", capture.output(print(st))))
    invisible(st)
}

### using both retstableR() and retstable()
set.seed(1)
alpha <- .2
V0 <- rgamma(2^12, 1)
set.seed(17); showSys.time(rET   <- retstable (alpha, V0)) ## method = default: here takes
							   ##	  983 times "MH",  17 x "LD"
set.seed(17); showSys.time(rET.H <- retstable (alpha, V0, method= "MH"))
set.seed(17); showSys.time(rET.D <- retstable (alpha, V0, method= "LD"))
set.seed(17); showSys.time(rET.R <- retstableR(alpha, V0))
T <- function(r) r^(1/8) # log() is too much
bp <- boxplot(T(rET), T(rET.H), T(rET.D), T(rET.R),
              notch=TRUE, col = "thistle")
(meds <- bp$stats[3,])

## "H0":  The 4 groups are not different -- here for the medians:
stopifnot(bp$conf[1,] < meds & meds < bp$conf[2,],
          bp$stats > 0,
          abs(bp$stats[2,] - 0.4035) < 0.007, ## first Quartiles
          abs(bp$stats[4,] - 0.8085) < 0.006) ## 3rd   Quartiles
