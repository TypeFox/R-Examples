
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
# FUNCTION:             PHASE SPACE REPRESENTATION:
#  mutualPlot            Creates mutual information plot
#  falsennPlot           Creates false nearest neigbours plot
# FUNCTION:             NON STATIONARITY:
#  recurrencePlot        Creates recurrence plot
#  separationPlot        Creates space-time separation plot
# FUNCTION:             LYAPUNOV EXPONENTS:
#  lyapunovPlot          Maximum Lyapunov plot              
################################################################################


test.mutualPlot = 
function()
{  
    # Mutual Information Index:
    par(mfrow = c(1, 1))
    lorentz = lorentzSim(
        times = seq(0, 40, by = 0.01), 
        parms = c(sigma = 16, r = 45.92, b = 4), 
        start = c(-14, -13, 47), 
        doplot = FALSE) 
    mutualPlot(x = lorentz[, 2], partitions = 16, lag.max = 20, doplot = TRUE) 
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.falsennPlot = 
function()
{  
    # False Nearest Neighbours:
    par(mfrow = c(1, 1))
    roessler = roesslerSim(
        times = seq(0, 100, by = 0.01), 
        parms = c(a = 0.2, b = 0.2, c = 8), 
        start = c(-1.894, -9.92, 0.025), 
        doplot = FALSE)
    falsennPlot(x = roessler[, 2], m = 6, d = 8, t = 180, eps = 1, rt = 3)
    abline(h = 0, col = "grey")
    grid()

   
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.recurrencePlot = 
function()
{  
    # Recurrence Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    lorentz = lorentzSim(
        times = seq(0, 40, by = 0.01), 
        parms = c(sigma = 16, r = 45.92, b = 4), 
        start = c(-14, -13, 47), 
        doplot = FALSE) 
    recurrencePlot(lorentz[, 2], m = 3, d = 2, end.time = 800, eps = 3, 
        nt = 5, pch = '.', cex = 2)
    recurrencePlot(lorentz[, 3], m = 3, d = 2, end.time = 800, eps = 3, 
        nt = 5, pch = '.', cex = 2)
    recurrencePlot(lorentz[, 4], m = 3, d = 2, end.time = 800, eps = 3, 
        nt = 5, pch = '.', cex = 2)
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.separationPlot = 
function()
{          
    # Separation Plot:
    par(mfrow = c(1, 1))
    roessler = roesslerSim(
        times = seq(0, 100, by = 0.01), 
        parms = c(a = 0.2, b = 0.2, c = 8), 
        start = c(-1.894, -9.92, 0.025), 
        doplot = FALSE)
    separationPlot(roessler[, 2], m = 3, d = 8, idt = 1, mdt = 250)    
   
    # Return Value:
    return()    
}


################################################################################


test.lyapunovPlot = 
function()
{          
    # Lyapunov Plot:
    NA
   
    # Return Value:
    return()    
}
   

################################################################################
    
