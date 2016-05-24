
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
#   1999 - 2006, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             DESCRIPTION:
#  dgh                   Returns density for generalized hyperbolic DF
#  pgh                   Returns probability for generalized hyperbolic DF
#  qgh                   Returns quantiles for generalized hyperbolic DF
#  rgh                   Returns random variates for generalized hyperbolic DF
# FUNCTION:             DESCRIPTION:
#  dhyp                  Returns density for hyperbolic DF
#  phyp                  Returns probability for hyperbolic DF
#  qhyp                  Returns quantiles for hyperbolic DF
#  rhyp                  Returns random variates for hyperbolic DF
#  hypMode               Computes the hyperbolic mode
# FUNCTION:             DESCRIPTION:
#  dnig                  Returns density for inverse Gaussian DF
#  pnig                  Returns probability for for inverse Gaussian DF
#  qnig                  Returns quantiles for for inverse Gaussian DF
#  rnig                  Returns random variates for inverse Gaussian DF
# FUNCTION:             DESCRIPTION:
#  hypSlider             Displays hyperbolic distribution function
#  nigSlider             Displays normal inverse Gausssian distribution function
################################################################################


test.gh =
function()
{
    # gh() Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    test = fBasics:::.distCheck("gh",
        alpha = 1.3, beta = 0.3, delta = 1.7, mu = 0.2, lambda = 0.8,
        n = 2000, robust = FALSE)
    print(test)
    checkTrue(mean(test) == 1)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.hyp =
function()
{
    # hyp() Distribution - Parameterization 1:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    test = fBasics:::.distCheck("hyp",
        alpha = 1.2, beta = 0.2, delta = 1.9, mu = 0.1, pm = "1",
        n = 1000, robust = FALSE)
    print(test)
    checkTrue(mean(test) == 1)

    # hyp() Distribution - Parameterization 2:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    test = fBasics:::.distCheck("hyp",
        alpha = 0.9, beta = -0.3, delta = 1.4, mu = -0.1, pm = "2",
        n = 1000, robust = FALSE)
    print(test)
    checkTrue(mean(test) == 1)

    # hyp() Distribution - Parameterization 3:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    fBasics:::.distCheck("hyp",
        alpha = 0.9, beta = -0.3, delta = 1.4, mu = -0.1, pm = "3",
        n = 1000, robust = FALSE)
    print(test)
    checkTrue(mean(test) == 1)

    # hyp() Distribution - Parameterization 4:
    if (FALSE) {
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        fBasics:::.distCheck("hyp",
            alpha = 1.6, beta = -0.3, delta = 1.4, mu = 0.1, pm = "4",
            n = 1000, robust = FALSE)                                    # CHECK
        print(test)
        checkTrue(mean(test) == 1)
    }

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.nig =
function()
{
    # nig() Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    test = fBasics:::.distCheck("nig",
        alpha = 2.1, beta = 0.1, delta = 1.5, mu = -0.1,
        n = 1000, robust = FALSE)
    print(test)
    checkTrue(mean(test) == 1)


    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.hypSlider =
function()
{
    # Arguments ?
    #   hypSlider()

    # Try:
    # hypSlider()
    NA

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.nigSlider =
function()
{
    # Arguments ?
    #   nigSlider

    # Try:
    # nigSlider()
    NA

    # Return Value:
    return()
}


################################################################################

