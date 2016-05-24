
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA.


################################################################################
# FUNCTION:                 TIME SERIES TESTS
#  lmTest                   Linear Modelling Test, select from:
#   bgTest                   Breusch-Godfrey Test
#   bpTest                   Breusch-Pagan Test
#   dwTest                   Durbin-Watson Test
#   gqTest                   Goldfeld-Quandt Test
#   harvTest                 Harvey-Collier Test
#   hmcTest                  Harrison-McCabe Test
#   rainTest                 Rainbow Test
#   resetTest                Ramsey's RESET Test
# REQUIRES:
# lmtest
################################################################################


################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: lmtest
#  Title: Testing Linear Regression Models
#  Version: 0.9-3
#  Date: $Date: 2003/02/19 15:54:30 $
#  Author: Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>,
#    Achim Zeileis <zeileis@ci.tuwien.ac.at>, David Mitchell
#  Maintainer: Achim Zeileis <zeileis@ci.tuwien.ac.at>
#  Description: A collection of tests, data sets and examples
#    for diagnostic checking in linear regression models.
#  Depends: R (>= 1.4.0)
#  License: GPL
################################################################################


lmTest <-
    function(formula,
    method = c("bg", "bp", "dw", "gq", "harv", "hmc", "rain", "reset"),
    data = list(), ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Load Library:
    #   Here we use a BUILTIN ...
    #   require(lmtest)

    # Settings:
    method <- match.arg(method)

    # DW:
    if (method == "dw")
        ans <- lmtest::dwtest(formula = formula, data = data, ...)

    # BP:
    if (method == "bp")
        ans <- lmtest::bptest(formula = formula, data = data, ...)

    # GQ:
    if (method == "gq")
        ans <- lmtest::gqtest(formula = formula, data = data, ...)

    # HMC:
    if (method == "hmc")
        ans <- lmtest::hmctest(formula = formula, data = data, ...)

    # HARV:
    if (method == "harv")
        ans <- lmtest::harvtest(formula = formula, data = data, ...)

    # RAIN:
    if (method == "rain")
        ans <- lmtest::raintest(formula = formula, data = data, ...)

    # RESET:
    if (method == "reset")
        ans <- lmtest::reset(formula = formula, data = data, ...)

    # BG:
    if (method == "bg")
        ans <- lmtest::bgtest(formula = formula, data = data, ...)

    # Return Result:
    ans
}


# ******************************************************************************


dwTest <-
    function(formula, alternative = c("greater", "two.sided", "less"),
    iterations = 15, exact = NULL, tol = 1.0e-10, data = list())
{
    lmtest::dwtest(formula, alternative, iterations, exact, tol, data)
}


# ------------------------------------------------------------------------------


bpTest <-
    function(formula, varformula = NULL, studentize = TRUE, data = list())
{
    lmtest::bptest(formula, varformula, studentize, data)
}


# ------------------------------------------------------------------------------


gqTest <-
    function(formula, point=0.5, order.by = NULL, data = list())
{
    lmtest::gqtest(formula, point, order.by, data)
}


# ------------------------------------------------------------------------------


hmcTest <-
    function(formula, point = 0.5, order.by = NULL, simulate.p = TRUE,
    nsim = 1000, plot = FALSE, data = list())
{
    lmtest::hmctest(formula, point, order.by, simulate.p, nsim, plot, data)
}


# ------------------------------------------------------------------------------


harvTest <-
    function(formula, order.by = NULL, data = list())
{
    lmtest::harvtest(formula, order.by, data)
}


# ------------------------------------------------------------------------------


rainTest =
    function(formula, fraction = 0.5, order.by = NULL, center = NULL,
    data = list())
{
    lmtest::raintest(formula, fraction, order.by, center, data)
}


# ------------------------------------------------------------------------------


resetTest <-
    function(formula, power = 2:3, type = c("fitted", "regressor", "princomp"),
    data = list())
{
    lmtest::reset(formula, power, type, data)
}


# ------------------------------------------------------------------------------


bgTest <-
    function(formula, order = 1, type = c("Chisq", "F"), data = list())
{
    lmtest::bgtest(formula, order, type, data)
}

################################################################################
