
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                   DESCRIPTION:
#  assetsRiskReturnPlot        Displays risk-return diagram of assets
#  assetsNIGShapeTrianglePlot  Displays NIG Shape Triangle
################################################################################


assetsRiskReturnPlot <-
    function(x, col = "steelblue",
    percentage = FALSE, scale = 252, labels = TRUE, add = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays risk-return giagram of assets

    # Arguments:
    #   x - a multivariate 'timeSeries' object
    
    # Example:
    #   x = 100 * as.timeSeries(data(LPP2005REC)) 
    #   assetsRiskReturnPlot(x)

    # FUNCTION:

    # Compute Return and Risk:
    if (percentage) index = 100 else index = 1

    # Compute Return and Risk:
    y = as.matrix(x)

    # Sample:
    Risk1 = index*sqrt(scale)* colStdevs(y)
    Return1 = index*scale*colMeans(y)

    # Huber(s):
    mu2 = mu3 = s2 = s3 = NULL
    for (i in 1:ncol(y)) {
        MeanSd2 = MASS::huber(y[, i])
        mu2 = c(mu2, MeanSd2$mu)
        s2 = c(s2, MeanSd2$s)
        # MeanSd3 = MASS::hubers(y[, i])
        # mu3 = c(mu3, MeanSd3$mu)
        # s3 = c(s3, MeanSd3$s)
    }
    Risk2 = index*sqrt(scale)*s2
    Return2 = index*scale*mu2
    # Risk3 = index*sqrt(scale)*s3
    # Return3 = index*scale*mu3

    # Colors:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)

    # Create Graph Frame:
    riskRange = range(c(Risk1, Risk2))
    riskRange[1] = 0
    riskRange[2] = riskRange[2] + 0.10*diff(riskRange)
    returnRange = range(c(Return1, Return2))
    returnRange[1] = returnRange[1] - 0.10*diff(returnRange)
    returnRange[2] = returnRange[2] + 0.10*diff(returnRange)

    if (labels) {
        plot(x = riskRange, y = returnRange,
            xlab = "Risk", ylab = "Return", type = "n")
        mtext("Sample versus Robust Estimates", line = 0.5, cex = 0.7)
    } else {
        plot(x = riskRange, y = returnRange,
            xlab = "", ylab = "", type = "n")
    }

    # Add all Points:
    colNames = colnames(x)
    for (i in 1:length(Risk1)) {
        points(Risk1[i], Return1[i], col = col[i], cex = 1.5, ...)
        if (add) {
            points(Risk2[i], Return2[i], col = col[i], cex = 1.1, ...)
        }
        text(
            Risk1[i] + diff(riskRange/50),
            Return1[i] + diff(returnRange/50),
            colNames[i], adj = 0, col = col[i])
    }
    if (labels) grid(col = "darkgrey")

    # Result:
    result = rbind(Risk1, Risk2, Return1, Return2)

    # Return Value:
    invisible(result)
}


# ------------------------------------------------------------------------------


assetsNIGShapeTrianglePlot <-
    function(x, labels = TRUE, col = "steelblue", ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays NIG Shape Triangle

    # Arguments:
    #   x - a multivariate 'timeSeries' object
    
    # Example:
    #   x = 100 * as.timeSeries(data(LPP2005REC)) 
    #   assetsNIGShapeTrianglePlot(x)

    # FUNCTION:

    # Settings:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)
    colNames = colnames(x)

    # Shape Triangle:
    for (i in 1:n) {
        fit = nigFit(100*x[, i], doplot = FALSE, trace = FALSE)
        nigShapeTriangle(fit, add = as.logical(i-1), labels = labels,
            col = col[i], ...)
        par = fit@fit$estimate
        alpha = par[1]
        beta = par[2]
        delta = par[3]
        mu = par[4]
        zeta = 1/sqrt(1 + delta * sqrt(alpha^2 - beta^2))
        chi = zeta * (beta/alpha)
        text(chi+0.01, zeta-0.01, colNames[i], adj = 0, col = col[i])
    }

    # Return Value:
    invisible()
}


################################################################################

