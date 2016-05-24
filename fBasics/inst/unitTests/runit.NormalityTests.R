
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
# FUNCTION:             NORMALITY TESTS:
#  normalTest            Test suite for normality tests
#  ksnormTest            One sample Kolmogorov-Smirnov normality test
#  shapiroTest           Shapiro-Wilk normality test
#  jarqueberaTest        Jarque-Bera normality test
#  dagoTest              D'Agostino normality test
# FUNCTION:             FROM NORTEST PACKAGE:
#  adTest                Anderson-Darling normality test
#  cvmTest               Cramer-von Mises normality test
#  lillieTest            Lilliefors (Kolmogorov-Smirnov) normality test
#  pchiTest              Pearson chi-square normality test
#  sfTest                Shapiro-Francia normality test
# FUNCTION:             MORE TESTS ...
#  runsTest              Runs test for detecting non-randomness [tseries]
#  gofnorm               Reports on several tests of normality
# FUNCTION ADDON:       DESCRIPTION:
#  jbTable               Table of finite sample p values for the JB test
#  pjb                   Computes probabilities for the Jarque Bera Test
#  qjb                   Computes quantiles for the Jarque Bera Test
#  jbTest                Performs finite sample adjusted JB LM and ALM test
################################################################################


test.NormalityTests =
function()
{
    # Normal Data:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    X = rnorm(50)

    TEST = ksnormTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = shapiroTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = jarqueberaTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    if(try(require(akima))) {
        TEST = jbTest(X)
        print(TEST)
        checkIdentical(as.character(class(TEST)), "fHTEST") }

    TEST = dagoTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = adTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = cvmTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = lillieTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = pchiTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = sfTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.NormalityTests.MSFT =
function()
{
    X = returnSeries(MSFT)[, 1]

    TEST = ksnormTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = shapiroTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = jarqueberaTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")


    if(try(require(akima))) {
        TEST = jbTest(X)
        print(TEST)
        checkIdentical(as.character(class(TEST)), "fHTEST")}

    TEST = dagoTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = adTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = cvmTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = lillieTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = pchiTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    TEST = sfTest(X)
    print(TEST)
    checkIdentical(as.character(class(TEST)), "fHTEST")

    # Return Value:
    return()
}


################################################################################


