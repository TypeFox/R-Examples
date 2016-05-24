
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


################################################################################
# Function:             DESCRIPTION:
# .varpiePlot            Plots a pie with variable radii
# .varpieDemo            Creates some demo plots
################################################################################


.varpiePlot <- 
function(x, radius = 0.8, labels = names(x), edges = 200, 
    clockwise = FALSE, init.angle = if (clockwise) 90 else 0, 
    density = NULL, angle = 45, col = NULL, border = NULL, 
    lty = NULL, main = NULL, ...) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots a pie with variable radii.
    
    # Note:
    #   The idea is given in:
    #   http://addictedtor.free.fr/packages/A2R/lastVersion/R/spieChart.R
  
    # Source:
    #   Package: A2R
    #   Type: Package
    #   Title: Romain Francois misc functions
    #   Version: 0.0-4
    #   Date: 2005-08-05
    #   Author: Romain Francois <francoisromain@free.fr>
    #   Maintainer: Romain Francois <francoisromain@free.fr>
    #   Description: Some functions of Romain Francois collection
    #   License: GPL version 2 or newer
    #   URL: http://addictedtor.free.fr/packages
  
    # FUNCTION:
    
    # Pie:
    if (length(radius) == 1) {
        # Normal Pie Plot:
        pie(x =x, labels = labels, edges = edges, radius = radius, 
            clockwise = clockwise, init.angle = init.angle, 
            density = density, angle = angle, col = col, 
            border = border, lty = lty, main = main, ...) 
    } else {
        # Variable Pie Plot:
        radius = radius/max(radius)
        X = x
        if (!is.numeric(x) || any(is.na(x) | x < 0)) 
            stop("'x' values must be positive.")
        if (is.null(labels)) {
            labels = as.character(1:length(x))
        } else {
            labels = as.graphicsAnnot(labels)
        }
        x = c(0, cumsum(x)/sum(x))
        dx = diff(x)
        nx = length(dx)
        plot.new()
        pin = par("pin")
        xlim = ylim = c(-1, 1)
        if (pin[1] > pin[2]) {
            xlim = (pin[1]/pin[2]) * xlim
        } else {
            ylim = (pin[2]/pin[1]) * ylim
        }
        plot.window(xlim, ylim, "", asp = 1)
        if (is.null(col)) col = rainbow(n)
        
        if (is.null(border)) border = "darkgrey"
        col = if (is.null(density)) col else par("fg")
        pie(x = X, col = "white", lty = NULL, 
            radius = max(radius), border = border, labels = labels, ...)
        phi = seq(0, 2*pi, length = 361)
        for (i in 1:length(radius))
            lines(radius[i]*cos(phi), radius[i]*sin(phi), col = border)    
                    
        col = rep(col, length.out = nx)
        border = rep(border, length.out = nx)
        lty = rep(lty, length.out = nx)
        angle = rep(angle, length.out = nx)
        density = rep(density, length.out = nx)
        twopi = if (clockwise) -2 * pi else 2 * pi
        t2xy = function(t, i) {
            t2p = twopi * t + init.angle * pi/180
            list(x = radius[i] * cos(t2p), y = radius[i] * sin(t2p))
        }
        
        for (i in 1:nx) {
            n = max(2, floor(edges * dx[i]))
            P = t2xy(seq(x[i], x[i + 1], length = n), i)
            polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
                border = border[i], col = col[i], lty = lty[i])
            P = t2xy(mean(x[i + 0:1]), i)
            # lab = as.character(labels[i])
            # if (!is.na(lab) && nchar(lab)) {
            #    lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y, col = "red")
            #    text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
            #       col = "red", adj = ifelse(P$x < 0, 1, 0), ...) }
        }
        title(main = main, ...)
        
        
    }
    
    # Return Values:
    invisible(NULL)
}


# ------------------------------------------------------------------------------


.varpieDemo <- 
function()
{
    # Description:
    #    Creates some demo plots
    
    # FUNCTION:
    
    # Demos:
    par (ask = TRUE)
     
    n = 64; phi = seq(0, pi/4, length = n)
    X = Y = runif(32)
    .varpiePlot(X, radius = Y, labels = "", border = "white")
    points(0, 0, pch = 19, cex = 4)
    
    
    n = 32; phi = seq(0, pi/4, length = n)
    X = abs(cos(phi)); Y = abs(sin(phi))
    .varpiePlot(X, radius = Y, labels = "", border = "white")
    points(0, 0, pch = 19, cex = 4)
    
    
    n = 64; phi = seq(0, pi, length = n)
    X = abs(cos(phi)); Y = abs(sin(phi))
    .varpiePlot(X, radius = Y, labels = "", border = "white", init.angle = -90)
    points(0, 0, pch = 19, cex = 4)
    
    
    n = 64; phi = seq(0, pi, length = n)
    X = abs(cos(phi)); Y = abs(sin(phi))
    .varpiePlot(X, radius = Y, labels = "", col = heat.colors(n), 
        border = "white", init.angle = -90)
    points(0, 0, pch = 19, cex = 4)
    
}   


################################################################################

