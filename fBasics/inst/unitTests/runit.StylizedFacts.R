
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
#  acfPlot               Displays autocorrelations function plot
#  pacfPlot              Displays partial autocorrelation function plot
#  teffectPlot           Estimates and plots the Taylor effect
#  lmacfPlot             Estimates and plots the long memory ACF
#  lacfPlot              Plots lagged autocorrelations
#  .logpdfPlot           Returns a pdf plot on logarithmic scale(s)
#  .qqgaussPlot          Returns a Gaussian quantile-quantile plot
#  scalinglawPlot        Evaluates and plots scaling law behavior
################################################################################


test.acfPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # acfPlot -
    acfPlot(x = msft.ret)

    # acfPlot -
    acfPlot(x = msft.ret, labels = FALSE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.pacfPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # pacfPlot -
    pacfPlot(x = msft.ret)

    # pacfPlot -
    pacfPlot(x = msft.ret, labels = FALSE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.teffectPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # teffectPlot -
    teffectPlot(x = msft.ret)

    # teffectPlot -
    teffectPlot(x = msft.ret, labels = FALSE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.lmacfPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # lmacfPlot -
    ## lmacfPlot(x = abs(msft.ret), type = "acf")
    ## lmacfPlot(x = abs(msft.ret), type = "hurst")
    # ... CHECK ACF OF RETURNS

    # lmacfPlot -
    ## lacfPlot(x = msft, n = 4, type = "values")                   ## CHECK !!!

    # lmacfPlot -
    ## lmacfPlot(x = abs(msft.ret), type = "acf", labels = FALSE)
    ## lmacfPlot(x = abs(msft.ret), type = "hurst", labels = FALSE)
    # ... CHECK ACF OF RETURNS

    # lmacfPlot -
    ## lacfPlot(x = msft, n = 4, labels = FALSE, type = "values")   ## CHECK !!!

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.lacfPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.logpdfPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # logpdfPlot -
    fBasics:::.logpdfPlot(x = msft.ret, labels = FALSE)
    fBasics:::.logpdfPlot(x = msft.ret, type = "log-log")
    # ... CHECK WARNINGS
    # ... CHECK COLORS

    # logpdfPlot -
    fBasics:::.logpdfPlot(x = msft.ret, labels = FALSE)
    fBasics:::.logpdfPlot(x = msft.ret, type = "log-log", labels = FALSE)
    # ... CHECK WARNINGS
    # ... CHECK COLORS

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.qqgausPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # qqgaussPlot -
    fBasics:::.qqgaussPlot(x = msft.ret)

    # qqgaussPlot -
    fBasics:::.qqgaussPlot(x = msft.ret, labels = FALSE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.scalinglawPlot =
function()
{
    # MSFT Data:
    msft.dat = MSFT
    msft = msft.dat[, 1]
    msft.vol = msft.dat[ , 5]/10^6
    msft.ret = returns(msft)

    # Graph Frame:
    par(mfrow = c(1, 1))

    # scalinglawPlot -
    scalinglawPlot(x = msft.ret, span = 4)
    # ... CHECK COLORS

    # scalinglawPlot -
    scalinglawPlot(x = msft.ret, span = 4, labels = FALSE)
    # ... CHECK COLORS

    # Return Value:
    return()
}


################################################################################

