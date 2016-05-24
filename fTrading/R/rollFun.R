
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  rollFun                   Compute Rolling Function Value
#   rollMean                  Compute Rolling Mean
#   rollVar                   Compute Rolling Variance
#   rollMin                   Compute Rolling Minimum
#   rollMax                   Compute Rolling Maximum
################################################################################


rollFun =
function(x, n, trim = TRUE, na.rm = FALSE, FUN, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute rolling function value

    # Arguments:
    #   x - an univariate "timeSeries" object or a numeric vector.
    #   n - an integer specifying the number of periods or
    #       terms to use in each rolling/moving sample.
    #   trim - a logical flag: if TRUE, the first n-1 missing values in
    #       the returned object will be removed; if FALSE, they will
    #       be saved in the returned object. The default is TRUE.
    #   na.rm - a logical flag: if TRUE, missing values in x will be
    #       removed before computation. The default is FALSE.
    #   FUN - the rolling function, arguments to this function can be
    #       passed through the \code{\dots} argument.

    # FUNCTION:

    # Transform:
    x.orig = x
    if (is.timeSeries(x)) {
        stopifnot(isUnivariate(x))
        TS = TRUE
    } else {
        TS = FALSE
    }
    if (TS) {
        positions = x.orig@positions
        x = series(x.orig)[, 1]
    } else {
        x = as.vector(x.orig)
        names(x) = NULL
    }

    # Remove NAs:
    if (na.rm) {
        if (TS) positions = positions[!is.na(x)]
        x = as.vector(na.omit(x))
    }

    # Roll FUN:
    start = 1
    end = length(x)-n+1
    m = x[start:end]
    if (n > 1) {
        for (i in 2:n) {
            start = start + 1
            end = end + 1
            m = cbind(m, x[start:end])
        }
    } else {
        m = matrix(m)
    }

    # Result:
    ans = apply(m, MARGIN = 1, FUN = FUN, ...)

    # Trim:
    if (!trim)
        ans = c(rep(NA, (n-1)), ans)
    if (trim & TS)
        positions = positions[-(1:(n-1))]

    # Back to timeSeries:
    if (TS) {
        ans = timeSeries(as.matrix(ans), positions, recordIDs = data.frame(),
            units = x.orig@units, FinCenter = x.orig@FinCenter)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


## rollMean =
## function(x, n = 9, trim = TRUE, na.rm = FALSE)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Compute rolling mean

##     # Examples:
##     #
##     #   x = timeSeries(as.matrix(cumsum(rnorm(12))), timeCalendar(),
##     #       units = "rnorm",FinCenter = "GMT")
##     #   rollMean(x, n = 4, trim = FALSE, na.rm = FALSE)
##     #   rollMean(x, n = 4, trim = TRUE, na.rm = FALSE)
##     #
##     #   series(x)[8, ] = NA
##     #   rollMean(x, n = 4, trim = FALSE, na.rm = FALSE)
##     #   rollMean(x, n = 4, trim = FALSE, na.rm = TRUE)
##     #   rollMean(x, n = 4, trim = TRUE, na.rm = TRUE)

##     # FUNCTION:

##     # Roll Mean:
##     rmean = rollFun(x = x, n = n, trim = trim, na.rm = na.rm, FUN = mean)

##     # Return Value:
##     rmean
## }


# ------------------------------------------------------------------------------


rollVar  =
function(x, n = 9, trim = TRUE, unbiased = TRUE, na.rm = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute rolling variance

    # FUNCTION:

    # Handle Time Series:
    if (is.timeSeries(x)) TS = TRUE else TS = FALSE

    # Roll Var:
    rvar = rollFun(x = x, n = n, trim = trim, na.rm = na.rm, FUN = var)

    # Unbiased ?
    if (!unbiased) {
        if (TS) {
            series(rvar) = (series(rvar) * (n-1))/n
        } else {
            rvar = (rvar * (n-1))/n
        }
    }

    # Return Value:
    rvar
}


# ------------------------------------------------------------------------------


## rollMax  =
## function(x, n = 9, trim = TRUE, na.rm = FALSE)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Compute rolling maximum

##     # FUNCTION:

##     # Roll Max:
##     rmax = rollFun(x = x, n = n, trim = trim, na.rm = na.rm,  FUN = max)

##     # Return Value:
##     rmax
## }


## # ------------------------------------------------------------------------------


## rollMin  =
## function(x, n = 9, trim = TRUE, na.rm = FALSE)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Compute rolling function minimum

##     # FUNCTION:

##     # Roll Min:
##     rmin = rollFun(x = x, n = n, trim = trim, na.rm = na.rm,  FUN = min)

##     # Return Value:
##     rmin
## }


################################################################################
