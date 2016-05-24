
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
# FUNCTION:                METRICS:
#  .riskMetricsPlot  
#  .garch11MetricsPlot     
################################################################################


.riskMetricsPlot <-
    function(x, labels = TRUE, lambda = 0.94, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #   x - an univariate timesSeries object
    
    # FUNCTION:
    
    # Check:
    stopifnot(isUnivariate(x))
    
    # Units:
    units = colnames(x)

    # Filter:
    riskMetrics = sqrt(emaTA(x^2, 1-lambda))
    
    # Plot:
    if (labels) {
        plot(riskMetrics, type = "l", col = "steelblue",
            main = paste(units, "RiskMetrics[TM]"), 
            xlab = "Time", ylab = "Volatility", ...)
        abline(h = sd(riskMetrics), col = "grey")
        SD = paste("StDev =", round(sd(x), 3))
        mtext(text = SD, side = 4, adj = 0, col = "grey", cex = 0.7)
        grid()
    } else {
        plot(riskMetrics, main = "", xlab = "", ylab = "", ...)
    }
    
    # Return Value:
    invisible(riskMetrics)     
} 


# ------------------------------------------------------------------------------


.garch11MetricsPlot <- 
    function(x, labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #   x - an univariate timesSeries object
    
    # FUNCTION:
    
    # Check:
    stopifnot(isUnivariate(x))
    
    # Units:
    units = colnames(x)

    # Filter:
    fit = garchFit(~garch(1,1), x, trace = FALSE)
    garch11 = volatility(fit)
    attr(garch11, "fit") <- fit
    
    # Plot:
    if (labels) {
        plot(garch11, type = "l", col = "steelblue",
            main = paste(units, "GARCH11 Volatility"), 
            xlab = "Time", ylab = "Volatility", ...)
        abline(h = sd(x), col = "grey")
        SD = paste("StDev =", round(sd(x), 3))
        mtext(text = SD, side = 4, adj = 0, col = "grey", cex = 0.7)
        grid()
    } else {
        plot(garch11, main = "", xlab = "", ylab = "", ...)
    }
    
    # Return Value:
    invisible(garch11)    
}


################################################################################

