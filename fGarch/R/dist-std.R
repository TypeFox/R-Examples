
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


################################################################################
# FUNCTION:              DESCRIPTION:
#  dstd                   Density for the Student-t Distribution
#  pstd                   Probability function for the Student-t Distribution
#  qstd                   Quantile function for the Student-t Distribution
#  rstd                   Random Number Generator for the Student-t
################################################################################


dstd <-
function(x, mean = 0, sd = 1, nu = 5, log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute the density for the
    #   Student-t distribution.

    # FUNCTION:

    # Params:
    if (length(mean) == 3) {
        nu = mean[3]
        sd = mean[2]
        mean = mean[1]
    }  
    
    # Compute Density:
    s = sqrt(nu/(nu-2))
    z = (x - mean) / sd
    result = dt(x = z*s, df = nu) * s / sd

    # Log:
    if(log) result = log(result)
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


pstd <-
function (q, mean = 0, sd = 1, nu = 5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute the probability for the
    #   Student-t distribution.

    # FUNCTION:

    # Compute Probability:
    s = sqrt(nu/(nu-2))
    z = (q - mean) / sd
    result = pt(q = z*s, df = nu)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


qstd <-
function (p, mean = 0, sd = 1, nu = 5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute the quantiles for the
    #   Student-t distribution.

    # FUNCTION:

    # Compute Quantiles:
    s = sqrt(nu/(nu-2))
    result = qt(p = p, df = nu) * sd / s + mean

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


rstd <-
function(n, mean = 0, sd = 1, nu = 5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generate random deviates from the
    #   Student-t distribution.

    # FUNCTION:

    # Generate Random Deviates:
    s = sqrt(nu/(nu-2))
    result = rt(n = n, df = nu) * sd / s + mean

    # Return Value:
    result
}


################################################################################

