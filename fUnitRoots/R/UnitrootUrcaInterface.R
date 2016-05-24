
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
# FUNCTION:                URCA INTERFACE [PFAFF] UNIT ROOT TESTS:
#  urdfTest                 Performs Augmented Dickey-Fuller test for unit roots
#  urersTest                PerformsElliott-Rothenberg-Stock test for unit roots
#  urkpssTest               Performs KPSS unit root test for stationarity
#  urppTest                 Performs Phillips-Perron test for unit roots
#  urspTest                 Performs Schmidt-Phillips test for unit roots
#  urzaTest                 Performs Zivot-Andrews test for unit roots
################################################################################


# NOTE, THIS INTERFACE REQUIRES THE URCA PACKAGE !


urdfTest =
function(x, lags = 1, type = c("nc", "c", "ct"), doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Dickey-Fuller test for unit roots

    # Notes:
    #   Requires "urca" which is not part of this distribution
    #   Wraps -
    #       ur.df(y, type = c("none", "drift", "trend"), lags = 1)

    # FUNCTION:

    # Compute:
    x = as.vector(x)
    if (type[1] == "nc") Type = "none"
    if (type[1] == "c")  Type = "drift"
    if (type[1] == "ct") Type = "trend"

    urca = ur.df(x, type = Type, lags = lags)
    output = capture.output(summary(urca))[-(1:4)]
    for (i in 1:length(output)) output[i] = paste(" ", output[i])
    output = output[-length(output)][-3]

    # Test Results:
    ans = list(
        name = "ur.df",
        test = urca,
        output = output
    )

    # Plot:
    if (doplot) plot(urca)

    # Return Value:
    new("fHTEST",
        call = match.call(),
        data = list(x = x),
        test = ans,
        title = "Augmented Dickey-Fuller Unit Root Test",
        description = description()
    )
}


# ------------------------------------------------------------------------------


urersTest =
function(x, type = c("DF-GLS", "P-test"), model = c("constant", "trend"),
lag.max = 4, doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Elliott-Rothenberg-Stock test for unit roots

    # Notes:
    #   Requires "urca" which is not part of this distribution
    #   Wraps -
    #       ur.ers(y, type = c("DF-GLS", "P-test"),
    #       model = c("constant", "trend"), lag.max = 4)

    # FUNCTION:

    # Compute:
    x = as.vector(x)
    urca = ur.ers(x, type = type[1], model = model[1], lag.max = lag.max)
    output = capture.output(summary(urca))[-(1:4)]
    for (i in 1:length(output)) output[i] = paste(" ", output[i])
    output = output[-length(output)]
    if (type[1] == "DF-GLS") output = output[-(4:7)]

    # Test Results:
    ans = list(
        name = "ur.ers",
        test = urca,
        output = output
    )

    # Plot:
    if (doplot & type[1] == "DF-GLS") plot(urca)

    # Return Value:
    new("fHTEST",
        call = match.call(),
        data = list(x = x),
        test = ans,
        title = "Elliott-Rothenberg-Stock Unit Root Test",
        description = description()
    )
}


# ------------------------------------------------------------------------------


urkpssTest =
function(x, type = c("mu", "tau"), lags = c("short", "long", "nil"),
use.lag = NULL, doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   KPSS unit root test for stationarity

    # Notes:
    #   Requires "urca" which is not part of this distribution
    #   Wraps:
    #       ur.kpss(y, type = c("mu", "tau"), lags = c("short", "long", "nil"),
    #       use.lag = NULL)

    # FUNCTION:

    # Compute:
    x = as.vector(x)
    urca = ur.kpss(x, type = type[1], lags = lags[1], use.lag = use.lag)
    output = capture.output(summary(urca))[-(1:4)]
    output = output[-length(output)]
    for (i in 1:length(output)) output[i] = paste(" ", output[i])

    # Test Results:
    ans = list(
        name = "ur.kpss",
        test = urca,
        output = output
    )

    # Plot:
    if (doplot) plot(urca)

    # Return Value:
    new("fHTEST",
        call = match.call(),
        data = list(x = x),
        test = ans,
        title = "KPSS Unit Root Test",
        description = description()
    )
}


# ------------------------------------------------------------------------------


urppTest =
function(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"),
lags = c("short", "long"), use.lag = NULL, doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Phillips-Perron test for unit roots

    # Note:
    #   Requires "urca" which is not part of this distribution
    #   Wraps:
    #       ur.pp(x, type = c("Z-alpha", "Z-tau"), model = c("constant",
    #       "trend"), lags = c("short", "long"), use.lag = NULL)

    # FUNCTION:

    # Compute:
    x = as.vector(x)
    urca = ur.pp(x, type = type[1], model = model[1], lags = lags[1],
        use.lag = use.lag)
    output = capture.output(summary(urca))[-c(1:4, 7:10)]
    for (i in 1:length(output)) output[i] = paste(" ", output[i])
    output = output[-length(output)]

    # Test Results:
    ans = list(
        name = "ur.pp",
        test = urca,
        output = output
    )

    # Plot:
    if (doplot) plot(urca)

    # Return Value:
    new("fHTEST",
        call = match.call(),
        data = list(x = x),
        test = ans,
        title = "Phillips-Perron Unit Root Test",
        description = description()
    )
}


# ------------------------------------------------------------------------------


urspTest =
function(x, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4),
signif = c(0.01, 0.05, 0.10), doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Schmidt-Phillips test for unit roots

    # Note:
    #   Requires "urca" which is not part of this distribution
    #   Wraps:
    #       ur.sp(y, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4),
    #       signif = c(0.01, 0.05, 0.1))

    # FUNCTION:

    # Compute:
    x = as.vector(x)
    urca = ur.sp(x, type = type[1], pol.deg = pol.deg[1], signif = signif[1])
    output = capture.output(summary(urca))[-(1:8)]
    output = output[-length(output)]
    for (i in 1:length(output)) output[i] = paste(" ", output[i])

    # Test Results:
    ans = list(
        name = "ur.pp",
        test = urca,
        output = output
    )

    # Plot:
    if (doplot) plot(urca)

    # Return Value:
    new("fHTEST",
        call = match.call(),
        data = list(x = x),
        test = ans,
        title = "Schmidt-Phillips Unit Root Test",
        description = description()
    )
}


# ------------------------------------------------------------------------------


urzaTest =
function(x, model = c("intercept", "trend", "both"), lag = 2, doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Zivot-Andrews test for unit roots

    # Note:
    #   Requires "urca" which is not part of this distribution
    #   Wraps:
    #       ur.za(y, model = c("intercept", "trend", "both"), lag)

    # FUNCTION:

    # Compute:
    x = as.vector(x)
    urca = ur.za(x, model = model[1], lag = lag)
    output = capture.output(summary(urca))[-(1:8)]
    output = output[-length(output)]
    for (i in 1:length(output)) output[i] = paste(" ", output[i])

    # Test Results:
    ans = list(
        name = "ur.pp",
        test = urca,
        output = output
    )

    # Plot:
    if (doplot) plot(urca)

    # Return Value:
    new("fHTEST",
        call = match.call(),
        data = list(x = x),
        test = ans,
        title = "Zivot & Andrews Unit Root Test",
        description = description()
    )
}


################################################################################
