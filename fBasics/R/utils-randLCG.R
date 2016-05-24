
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
# FUNCTION:                 PORTABLE INNOVATIONS:
#  set.lcgseed               Sets initial random seed
#  get.lcgseed               Gets the current valus of the random seed
# FUNCTION:                 PORTABLE RVS GENERATORS:
#  runif.lcg                 Generates uniform linear congruational rvs
#  rnorm.lcg                 Generates normal linear congruational rvs
#  rt.lcg                    Generates Student-t linear congruational rvs
################################################################################


################################################################################
#  set.lcgseed               Sets initial random seed
#  get.lcgseed               Gets the current valus of the random seed


# .lcg.seed = 4711


# ------------------------------------------------------------------------------


set.lcgseed <-
function(seed = 4711)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets the random seed for the linear congruential
    #   random number generator

    # Notes:
    #   A Simple Portable Random Number Generator for Use in R and Splus
    #   Use this generator only for comparisons of Programs in R and Splus !!!
    #   Method: A linear congruential generator with
    #       LCG(a=13445, c=0, m=2^31-1, X0)
    #       Note, this is a random number generator which passes the bitwise
    #       randomness test.
    #   Reference:
    #       http://csep1.phy.ornl.gov/rn/node13.html
    #       N. S. Altman. ``Bitwise Behavior of Random Number Generators,''
    #       SIAM J. Sci. Stat. Comput., 9(5), September, pps. 941-949, 1988

    #   Example:
    #       set.lcgseed(4711)
    #       cbind(runif.lcg(100), rnorm.lcg(100), rt.lcg(100, df=4))

    # FUNCTION:

    # Return Value:
    # .lcg.seed <<- seed
    setRmetricsOptions(lcg.seed = seed)
}


# ------------------------------------------------------------------------------


get.lcgseed <-
function()
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the random seed for the linear congruential
    #   random number generator

    # FUNCTION

    # Return Value:
    getRmetricsOptions("lcg.seed")
}


################################################################################
#  runif.lcg                 Uniform linear congruational generator
#  rnorm.lcg                 Normal linear congruational generator
#  rt.lcg                    Student-t linear congruational generator


runif.lcg <-
function(n, min = 0, max = 1)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #    A linear congruential generator for uniform distributed
    #    random numbers

    # Notes:
    #    Important - Use this generator only for comparisons of
    #      Programs in R and SPlus (and not for production) !!!
    #    Portable Random Numbers:
    #      A linear congruential generator
    #      LCG(a=13445, c=0, m=2^31-1, X0)
    #      This is a random number generator which
    #      passes the bitwise randomness test

    # References:
    #    http://csep1.phy.ornl.gov/rn/node13.html
    #    N. S. Altman. ``Bitwise Behavior of Random Number Generators,''
    #    SIAM J. Sci. Stat. Comput., 9(5), September, pps. 941-949, 1988

    # FUNCTION:

    # Initialize:
    # if(!exists(".lcg.seed")) .lcg.seed <<- 4711

    # Generate:
    r.lcg = rep(0, times = n)
    a = 13445
    c = 0
    m = 2^31-1
    for (i in 1:n) {
        lcg.seed <- getRmetricsOptions("lcg.seed")
        setRmetricsOptions(lcg.seed = (a * lcg.seed + c) %% m)
        r.lcg[i] = getRmetricsOptions("lcg.seed") / m }
    r.lcg = (max-min)*r.lcg + min

    # Return Value:
    r.lcg
}


# ------------------------------------------------------------------------------


rnorm.lcg <-
function(n, mean = 0, sd = 1)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #    A linear congruential generator for normal distributed
    #    random numbers

    # FUNCTION:

    # This is slow, but portable between R and SPlus
    (qnorm(runif.lcg(n = n, min = 0, max = 1)) - mean)/sd^2
}


# ------------------------------------------------------------------------------


rt.lcg <-
function(n, df)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #    A linear congruential generator for Sudent-t distributed
    #    random numbers

    # FUNCTION:

    # This is slow, but portable between R and SPlus
    qt(runif.lcg(n = n, min = 0, max = 1), df = df)
}


################################################################################

