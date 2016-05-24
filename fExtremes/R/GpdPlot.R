
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
# METHODS:                PRINT, PLOT, AND SUMMARY:
#  plot.fGPDFIT            S3 Plot Method for object of class "fGPDFIT"
#  .gpd1Plot                Empirical Distribution Plot
#  .gpd2Plot                Tail of Underlying Distribution
#  .gpd3Plot                Scatterplot of GPD Residuals
#  .gpd4Plot                Quantile-Quantile Plot of GPD Residuals
################################################################################


plot.fGPDFIT =
function(x, which = "ask", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for objects of class 'fGPDFIT'

    # Example:
    #   x = as.timeSeries(danishClaims); plot(gpdFit(x, 4), "ask")
    
    # FUNCTION:
            
    # Plot:
    interactivePlot(
        x = x,
        choices = c(
            "Excess Distribution",
            "Tail of Underlying Distribution",
            "Scatterplot of Residuals", 
            "QQ-Plot of Residuals"),
        plotFUN = c(
            ".gpd1Plot", 
            ".gpd2Plot", 
            ".gpd3Plot",
            ".gpd4Plot"),
        which = which)
                    
    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.gpd1Plot =
function(x, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Empirical Distribution Plot
    
    # Arguments:
    #   x - an object of class fGPDFIT
    #   labels - a logical flag. Should labels be printed?
    
    # FUNCTION:
    
    # Data:
    extend = 1.5
    u = x@parameter$u
    data = as.vector(x@data$exceedances)
    sorted = sort(data)
    shape = xi = x@fit$par.ests["xi"]
    scale = beta = x@fit$par.est["beta"]
    ypoints = ppoints(sorted)
    U = max(sorted)*extend
    z = qgpd(seq(0, 1, length = 1000), xi, u, beta)
    z = pmax(pmin(z, U), u)
    y = pgpd(z, xi, u, beta) 
    
    # Labels:
    if (labels) {
        xlab = "Fu(x-u)"
        ylab = "x [log Scale]"
        main = "Excess Distribution"
    } else {
        xlab = ylab = main = ""
    }
    
    # Plot:
    plot(x = sorted, y = ypoints, 
        xlim = range(u, U), ylim = range(ypoints, y, na.rm = TRUE), 
        main = main, xlab = xlab, ylab = ylab, 
        log = "x", axes = TRUE, 
        col = "steelblue", pch = 19, ...)
    lines(z[y >= 0], y[y >= 0], col = "brown") 
    
    # Addon:
    if (labels) {
        u = signif (x@parameter$u, 3)
        text = paste("u =", u)
        mtext(text, side = 4, adj = 0, cex = 0.7)
        grid()
    }
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.gpd2Plot =
function(x, labels = TRUE, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Tail of Underlying Distribution
    
    # Arguments:
    #   x - an object of class fGPDFIT
    #   labels - a logical flag. Should labels be printed?
    
    # FUNCTION:
    
    # Settings:
    extend = 1.5
    u = x@parameter$u
    data = as.vector(x@data$x)
    sorted = sort(data[data > u])
    prob = x@fit$prob
    shape = xi = x@fit$par.ests["xi"]
    beta = x@fit$par.ests["beta"]
    scale = beta * (1-prob)^xi
    location = u - (scale*((1 - prob)^(-xi)-1))/xi

    # Labels:
    if (labels) {
        xlab = "x [log scale]"
        ylab = "1-F(x) [log scale]"
        main = "Tail of Underlying Distribution"
    } else {
        xlab = ylab = main = ""
    }
    
    # Plot:
    U = max(data) * extend
    ypoints = ppoints(sorted)
    ypoints = (1 - prob) * (1 - ypoints)
    z = qgpd(seq(0, 1, length = 1000), xi, u, beta)
    z = pmax(pmin(z, U), u)
    y = pgpd(z, xi, u, beta)   
    y = (1 - prob) * (1 - y)
    plot(x = sorted, y = ypoints, 
        xlim = range(u, U), ylim = range(ypoints, y[y>0], na.rm = TRUE), 
        main = main, xlab = xlab, ylab = ylab, 
        log = "xy", axes = TRUE, 
        col = "steelblue", pch = 19, ...)
        
    # Line:
    lines(z[y >= 0], y[y >= 0], col = "brown")    
    if (labels) grid()
    
    # Return Value:
    invisible(list(x = sorted, y = ypoints))
}


# ------------------------------------------------------------------------------


.gpd3Plot =
function(x, labels = TRUE, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Scatterplot of GPD Residuals
    
    # Arguments:
    #   x - an object of class fGPDFIT
    #   labels - a logical flag. Should labels be printed?
    
    # FUNCTION:
    
    # Residuals:
    residuals = x@residuals
    
    # Labels:
    if (labels) {
        ylab = "Residuals"
        xlab = "Ordering"
        main = "Scatterplot of Residuals"
    } else {
        xlab = ylab = main = ""
    }
    
    # Plot:
    plot(residuals, 
       main = main, ylab = ylab, xlab = xlab, 
       col = "steelblue", pch = 19, ...)
    lines(lowess(1:length(residuals), residuals), col = "brown") 
    if (labels) grid()
    
    # Return Value:
    invisible(list(x = 1:(length(residuals)), y = residuals))
}


# ------------------------------------------------------------------------------


.gpd4Plot = 
function(x, labels = TRUE, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Quantile-Quantile Plot of GPD Residuals
    
    # Arguments:
    #   x - an object of class fGPDFIT
    #   labels - a logical flag. Should labels be printed?
    
    # FUNCTION:
    
    # Data:
    data = x@residuals
    sorted = sort(data)
    
    # Labels:
    if (labels) {
        xlab = "Ordered Data"
        ylab = "Exponential Quantiles"
        main = "QQ-Plot of Residuals"
    } else {
        xlab = ylab = main = ""
    }
    
    # Plot:
    y = qexp(ppoints(data))
    plot(x = sorted, y = y, 
        main = main, xlab = xlab, ylab = ylab, 
        col = "steelblue", pch = 19, ...)
    abline(lsfit(sorted, y), col = "brown")
    if (labels) grid()
    
    # Return Value:
    invisible(list(x = sorted, y = y))
}


################################################################################

