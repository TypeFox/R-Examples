
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


test.mutualPlotGallery = 
function()
{  
    # Mutual Information Index:
    lorentz = lorentzSim(
        times = seq(0, 40, by = 0.01), 
        parms = c(sigma = 16, r = 45.92, b = 4), 
        start = c(-14, -13, 47), 
        doplot = FALSE)
        
    # Plot: 
    par(mfrow = c(1, 1))
    mutualPlot(x = lorentz[, 2], partitions = 16, lag.max = 20, doplot = TRUE) 
    mtext("Lorentz Map", line = 0.5, cex = 0.7)
    mtext(paste("times=seq(0,40,by=0.01) | parms=c(sigma=16,r=45.92,b=4) |",
        "start=c(-14,-13,47)"), side = 4, adj = 0, col = "darkgrey", cex = 0.7)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.falsennPlotGallery = 
function()
{  
    # False Nearest Neighbours:
    roessler = roesslerSim(
        times = seq(0, 100, by = 0.01), 
        parms = c(a = 0.2, b = 0.2, c = 8), 
        start = c(-1.894, -9.92, 0.025), 
        doplot = FALSE)
    
    # Plot:
    par(mfrow = c(1, 1))
    falsennPlot(x = roessler[, 2], m = 6, d = 8, t = 180, eps = 1, rt = 3)
    abline(h = 0, col = "grey")
    grid()
    mtext("Roessler Map", line = 0.5, cex = 0.7)
    mtext(paste("times=seq(0,100,by=0.01) | parms=c(a=0.2, b=0.2, c=8) |",
        "start=c(-1.894,-9.92,0.025)"), side = 4, adj = 0, 
        col = "darkgrey", cex = 0.7)
       
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------