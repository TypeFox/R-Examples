
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
# FUNCTION:               ADDITIONAL PLOTS:
#  gpdTailPlot             Plots Tail Estimate From GPD Model
#  gpdQuantPlot            Plots of GPD Tail Estimate of a High Quantile
#  gpdShapePlot            Plots for GPD Shape Parameter
#  gpdQPlot                Adds Quantile Estimates to plot.gpd
#  gpdSfallPlot            Adds Expected Shortfall Estimates to a GPD Plot
#  gpdRiskMeasures         Calculates Quantiles and Expected Shortfalls
# FUNCTION:               NEW STYLE FUNCTIONS:
#  tailPlot                Plots GPD VaR and Expected Shortfall risk
#  tailSlider              Interactive view to find proper threshold value
#  tailRiskMeasures        Calculates VaR and Expected Shortfall risks
################################################################################


test.gpdTailPlot = 
function()
{
    # Artificial Data Set:
    x = gpdSim(seed = 1985)
    fit = gpdFit(x)
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    gpdTailPlot(fit)
    
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    fit = gpdFit(x)
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    gpdTailPlot(fit)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gpdQuantPlot = 
function()
{
    # Artificial Data Set:
    x = gpdSim(seed = 1985)
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    gpdQuantPlot(x)
    
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    fit = gpdFit(x)
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    gpdQuantPlot(x)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gpdShapePlot = 
function()
{
    # Artificial Data Set:
    x = gpdSim(seed = 1985)
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    gpdShapePlot(x)
    
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    gpdShapePlot(x)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gpdQPlot = 
function()
{
    # Artificial Data Set:
    x = gpdSim(seed = 1985)
    fit = gpdFit(x)
    tp = gpdTailPlot(fit)
    gpdQPlot(tp)
    
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    fit = gpdFit(x, u =10)
    tp = gpdTailPlot(fit)
    gpdQPlot(tp)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gpdSfallPlot = 
function()
{
    # Artificial Data Set:                                 
    x = gpdSim(seed = 1985)
    fit = gpdFit(x)
    ### tp = gpdTailPlot(fit)                                            # CHECK
    ### gpdSfallPlot(tp)                                                 # CHECK
    
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    fit = gpdFit(as.vector(x), u =10)
    ### tp = gpdTailPlot(fit)                                            # CHECK
    ### gpdSfallPlot(tp)                                                 # CHECK

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.tailPlot = 
function()
{
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    fit = gpdFit(x, u = 10)
    ### tailPlot(fit)                                                    # CHECK

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.tailSlider = 
function()
{
    # Danish Fire Claims:
    # x = as.timeSeries(data(danishClaims))
    # tailSlider(x)
    NA
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.tailRisk = 
function()
{
    # Danish Fire Claims:
    x = as.timeSeries(data(danishClaims))
    fit = gpdFit(x, u = 10)
    tailRisk(fit)

    # Return Value:
    return()    
}


################################################################################

