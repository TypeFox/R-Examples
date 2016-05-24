# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C) for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# and
# Copyright (C) 2010--2012 Martin Maechler, ETH Zurich
#
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTIONS:            STABLE DISTRIBUTION:
#  stableMode            Computes stable mode
#  dstable               Returns density for stable DF
#  pstable               Returns probabilities for stable DF
#  qstable               Returns quantiles for stable DF
#  rstable               Returns random variates for stable DF
################################################################################

## FIXME / TODO:  need unit tests  for   stableMode()  !!


if(do.stable.rUnitTest <- require("fBasics") &&
   ## need the newer distCheck():
   packageDescription("fBasics")$Version >= package_version("2110.79"))
{
    ## fBasics:  for  .distCheck()
    distCheck <- fBasics:::.distCheck
    environment(distCheck) <- asNamespace("stabledist")
    ## and re-attach "stabledist" as its contents is now masked by fBasics::dstable:
    if((P <- "package:stabledist") %in% search())
	detach(P, character.only=TRUE)
    stopifnot(require("stabledist"),
	      ## check that indeed we get stabledist's functions, not fBasics:
	      identical(dstable, stabledist::dstable))
}
source(system.file("test-tools-1.R", package = "Matrix"))
					#-> identical3(), showProc.time(),...
(doExtras <- stabledist:::doExtras())
n.check <- if(doExtras) 1000 else 64

test.stableS0 <- function(n = n.check)
{
    if (do.stable.rUnitTest) {
        ## "TODO" in distCheck() -- use  'tol = .005'  in newer versions
        # stable - Parameterization S0:
        test <- distCheck("stable", n=n, alpha = 1.8, beta = 0.3)
        print(test)
        ## the 3rd test -- matching (mean, var)  typically fails for stable -- as Var(.) == Inf !
        checkTrue(all(test[1:2]))

        # stable - Parameterization S0:
        test <- distCheck("stable", n=n, alpha = 1.2, beta = -0.3)
        print(test)
        checkTrue(all(test[1:2]))

      if(doExtras) {
        # stable - Parameterization S0:
        test <- distCheck("stable", n=n, alpha = 0.6, beta = 0)
        print(test)
        checkTrue(all(test[1:2]))
      }
    }

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.stableS1 <- function(n = n.check)
{
    if (do.stable.rUnitTest) {
      if(doExtras) {
        # stable - Parameterization S1:
        test <- distCheck("stable", n=n, alpha = 1.8, beta = 0.3, pm = 1)
        print(test)
        checkTrue(all(test[1:2]))
      }
        # stable - Parameterization S1:
        test <- distCheck("stable", n=n, alpha = 1.2, beta = -0.3, pm = 1)
        print(test)
        checkTrue(all(test[1:2]))

      if(doExtras) {
        # stable - Parameterization S1:
        test <- distCheck("stable", n=n, alpha = 0.6, beta = 0, pm = 1)
        print(test)
        checkTrue(all(test[1:2]))
      }
    }

    # Return Value:
    return()
}
##if(doExtras) test.stableS1 <- Tst.stableS1

# ------------------------------------------------------------------------------


Tst.stableS2 <- function(n = n.check)
{
    if (do.stable.rUnitTest) {
      if(doExtras) {
        # stable - Parameterization S2:
        test <- distCheck("stable", n=n, alpha = 1.8, beta = 0.3, pm = 2)
        print(test)
        checkTrue(all(test[1:2]))

        # stable - Parameterization S2:
        test <- distCheck("stable", n=n, alpha = 1.2, beta = -0.3, pm = 2)
        print(test)
        checkTrue(all(test[1:2]))
      }
        # stable - Parameterization S2:
        test <- distCheck("stable", n=n, alpha = 0.6, beta = 0, pm = 2)
        print(test)
        checkTrue(all(test[1:2]))
    }

    # Return Value:
    return()
}
if(doExtras) test.stableS2 <- Tst.stableS2

################################################################################

