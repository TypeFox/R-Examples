
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
# FUNCTION:             MDA ESTIMATORS:
#  hillPlot              Plot Hill's estimator
#  shaparmPlot           Pickands, Hill & Decker-Einmahl-deHaan Estimator
#   shaparmPickands      Auxiliary function called by shaparmPlot
#   shaparmHill           ... called by shaparmPlot
#   shaparmDehaan         ... called by shaparmPlot
################################################################################


test.hillPlot = 
function()
{
    #  hillPlot              Plot Hill's estimator
   
    # Graph Frame:
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    
    # Hill Plot:
    hillPlot(gevSim(n=1000), plottype = "alpha")
    hillPlot(gevSim(n=1000), plottype = "xi"); grid()
    
    # Don't Plot Return Value:
    hillPlot(gevSim(n=1000), plottype = "alpha", doplot = FALSE)
    hillPlot(gevSim(n=1000), plottype = "xi", doplot = FALSE); grid()   
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.shaparmPlot = 
function()
{
    #  shaparmPlot           Pickands, Hill & Decker-Einmahl-deHaan Estimator
   
    # Graph Frame:
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    
    # shaparmPlot(x, p = 0.01*(1:10), xiRange = NULL, alphaRange = NULL,
    #   doplot = TRUE, plottype = c("both", "upper"))

    # Graph Frame:
    par(mfcol = c(3, 2), cex = 0.7)
    par(ask = FALSE)
    shaparmPlot(as.timeSeries(data(bmwRet)))
    
    # Print (Results:
    shaparmPlot(as.timeSeries(data(bmwRet)), doplot = FALSE)
    
    # Tailored p:
    shaparmPlot(as.timeSeries(data(bmwRet)), p = 0.005*(2:20))
    
    
    # Return Value:
    return()    
}


################################################################################

