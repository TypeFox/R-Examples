
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


# ##############################################################################
# FUNCTION:               GPD SIMULATION:
#  gpdSim                  Simulates a GPD distributed process
# FUNCTION:               GPD PARAMETER ESTIMATION:
#  'fGPDFIT'               S4 class representation
#  gpdFit                  Fits Parameters of GPD distribution
# METHODS:                PRINT, PLOT, AND SUMMARY:
#  show.fGPDFIT            S4 Print Method for object of class "fGPDFIT"
#  plot.fGPDFIT            S3 Plot Method for object of class "fGPDFIT"
#  summary.fGPDFIT         S3 Summary Method for object of class "fGPDFIT"
################################################################################


test.gpdSim =
function()
{
    # Generate Artificial Data Set:
    x = gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000, seed = 4711)
    class(x)

    # Plot Series:
    par(mfrow = c(2, 1), cex = 0.7)
    par(ask = FALSE)
    seriesPlot(as.timeSeries(x))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.fGPDFIT =
function()
{
    # Slot names:
    slotNames("fGPDFIT")
    # [1] "call"        "method"      "parameter"   "data"        "fit"
    # [6] "residuals"   "title"       "description"

    # Return Value:
    return()
}



# ------------------------------------------------------------------------------


test.gpdFit =
function()
{
    # Generate Artificial Data Set:
    model = list(xi = -0.25, mu = 0, beta = 1)
    ts = gpdSim(model = model, n = 5000, seed = 4711)
    class(ts)

    # Transform As timeSeries:
    tS = as.timeSeries(ts)
    class(tS)

    # As numeric vector:
    x = as.vector(ts)
    class(x)

    # GPD Fit:
    # gpdFit(x, u = quantile(x, 0.95), type = c("mle", "pwm"),
    #   information = c("observed", "expected"), title = NULL,
    #   description = NULL, ...)

    # PWM Fit:
    fit = gpdFit(tS, u = min(series(tS)), "pwm")
    print(fit)
    fit = gpdFit(ts, u = min(ts), "pwm")
    print(fit)
    fit = gpdFit(x, u = min(x), "pwm")
    print(fit)

    # MLE Fit:
    fit = gpdFit(tS, u = min(series(tS)), "mle")
    print(fit)
    fit = gpdFit(ts, u = min(ts), "mle")
    print(fit)
    fit = gpdFit(x, u = min(x), "mle")
    print(fit)

    # Information:
    fit = gpdFit(tS, u = min(series(tS)), type = "mle", information = "observed")
    print(fit)
    fit = gpdFit(tS, u = min(series(tS)), type = "mle", information = "expected")
    print(fit)


    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.plot =
function()
{
    # Artificial Data Set:
    model = list(xi = -0.25, mu = 0, beta = 1)
    ts = gpdSim(model = model, n = 5000, seed = 4711)
    class(ts)

    # Fit:
    fit = gpdFit(ts, u = min(ts), type = "mle")
    print(fit)
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    plot(fit, which = "all")

    # Try:
    # plot(fit, which = "ask")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.summary =
function()
{
    # Artificial Data Set:
    model = list(xi = -0.25, mu = 0, beta = 1)
    ts = gpdSim(model = model, n = 5000, seed = 4711)
    class(ts)

    # Fit:
    fit = gpdFit(ts, u = min(ts), type = "mle")
    summary(fit, doplot = FALSE)

    # Return Value:
    return()
}


################################################################################

