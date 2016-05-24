
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
# FUNCTION:             DESCRIPTION:
#  dgld                 Returns density for Generalized Lambda DF
#  pgld                 Returns probability for Generalized Lambda DF
#  qgld                 Returns quantiles for Generalized Lambda DF
#  rgld                 Returns random variates for Generalized Lambda DF
################################################################################


dgld <-
function(x, lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8, log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns density for Generalized Lambda DF

    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter

    # Details:
    #   This a simple wrapper for the GLD using the RS parameterisation
    #   in Region 4. The parameter range is restricted to negative lambda3
    #   and lambda4 parameters

    # Note:
    #   Uses gld functions from package gld

    # Example:
    #   dgld( (-25:25)/10 )

    # FUNCTION:

    # Parameters:
    if (length(lambda1) == 4) {
       lambda4 = lambda1[4]
       lambda3 = lambda1[3]
       lambda2 = lambda1[2]
       lambda1 = lambda1[1]
    }

    # Check Parameters:
    stopifnot (lambda2 < 0)
    stopifnot (lambda3 < 0)
    stopifnot (lambda4 < 0)

    # Density:
    d = .dgld(x,
        lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, lambda4=lambda4,
        param="rs")

    # Log:
    if(log) d = log(d)

    # Return Value:
    d
}


# ------------------------------------------------------------------------------


pgld <-
function(q, lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns probability for Generalized Lambda DF

    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter

    # Note:
    #   Uses gld functions from package gld

    # Example:
    #   pgld( (-25:25)/10 )

    # FUNCTION:

    # Check Parameters:
    stopifnot (lambda2 < 0)
    stopifnot (lambda3 < 0)
    stopifnot (lambda4 < 0)

    # Probability:
    p = .pgld(q,
        lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, lambda4=lambda4,
        param="rs")

    # Return Value:
    p
}


# ------------------------------------------------------------------------------


qgld <-
function(p, lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles for Generalized Lambda DF

    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter

    # Note:
    #   Uses gld functions from package gld

    # Example:
    #   qgld((1:99)/100)

    # FUNCTION:

    # Check Parameters:
    stopifnot (lambda2 < 0)
    stopifnot (lambda3 < 0)
    stopifnot (lambda4 < 0)

    # Quantiles:
    q = .qgld(p,
        lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, lambda4=lambda4,
        param="rs")

    # Return Value:
    q
}


# ------------------------------------------------------------------------------


rgld <-
function(n = 100, lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random variates for Generalized Lambda DF

    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter

    # Note:
    #   Uses gld functions from package gld

    # FUNCTION:

    # Check Parameters:
    stopifnot (lambda2 < 0)
    stopifnot (lambda3 < 0)
    stopifnot (lambda4 < 0)

    # Random Variates:
    r = .rgld(n,
        lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, lambda4=lambda4,
        param="rs")

    # Return Value:
    r
}


################################################################################

