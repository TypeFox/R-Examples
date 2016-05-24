
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
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file 


################################################################################
# FUNCTION:            PLOT UTILITIES:
#  interactivePlot      Plots several graphs interactively
################################################################################


test.interactivePlot <- 
    function()
{
    # interactivePlot(x, choices = paste("Plot", 1:9), 
    #   plotFUN = paste("plot.", 1:9, sep = ""), which = "all", ...)
    
    # Test Plot Function:
    testPlot = function(x, which = "all", ...) {   
        # Plot Function and Addons:
        plot.1 <<- function(x, ...) plot(x, ...)      
        plot.2 <<- function(x, ...) acf(x, ...)
        plot.3 <<- function(x, ...) hist(x, ...)      
        plot.4 <<- function(x, ...) qqnorm(x, ...)
        # Plot:
        interactivePlot(x,
            choices = c("Series Plot", "ACF", "Histogram", "QQ Plot"),
            plotFUN = c("plot.1", "plot.2", "plot.3", "plot.4"),
            which = which, ...)       
        # Return Value:
        invisible()
    } 
    
    # Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    testPlot(rnorm(500))
    
    # Reorder 2:
    par(mfrow = c(2, 1), cex = 0.7)
    testPlot(rnorm(500), which = c(1, 2))
       
    # Try:
    # par(mfrow = c(1,1)); testPlot(rnorm(500), which = "ask")
        
    # Return Value:
    return()
}
       

################################################################################

