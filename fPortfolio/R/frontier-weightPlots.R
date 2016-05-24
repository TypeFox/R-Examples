
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                    DESCRIPTION:
#  .weightsWheel                Adds a pie of weights to frontier plot
#  .attributesWheel             Adds a pie of attributes to frontier plot
# FUNCTION:                    DESCRIPTION:
#  .notStackedWeightsPlot       Plots the not stacked weights of potfolio
#  .addlegend                   Adds legend to sliders
################################################################################


.weightsWheel <-
    function(object, piePos = NULL, pieR = NULL, pieOffset = NULL, ...)
{
    # A function implemented by Oliver Greshake and Diethelm Wuertz

    # Description:
    #   Adds a pie plot of weights for MV and CVaR Portfolios

    # Arguments:
    
    # Details:
    #   The default settings are:
    #   piePos - Position of tangency Portfolio
    #   pieR - 10% of the Risk Range: diff(range(targetRisk(object)))/10

    # FUNCTION:

    # Extract Coordinates:
    p = par()$usr/15
    dx = p[2]-p[1]
    dy = p[4]-p[3]

    # Pie Position:
    if(is.null(piePos)) {
        Data = getSeries(object)
        Spec = getSpec(object)
        Constraints = getConstraints(object)
        tg = getTargetReturn(tangencyPortfolio(Data, Spec, Constraints))
        ef = as.vector(getTargetReturn(object))
        piePos = which(diff(sign(ef-tg)) > 0)
    }

    # Pie Radius:
    if(is.null(pieR)) {
        pieR = c(1, 1)
    }

    # Pie Offset:
    if(is.null(pieOffset)) {
        pieOffset = c(-2*dx, 0)
    }

    # Plot Circle:
    weights = getWeights(object)[piePos, ]
    nWeights = length(weights)
    Sign = rep("+", nWeights)
    Sign[(1:nWeights)[weights < 0]] = "-"
    x = getTargetRisk(object)[piePos]
    y = getTargetReturn(object)[piePos]
    phi =  seq(0, 2*pi, length = 360)
    X = x + pieOffset[1] + pieR[1] * sin(phi) * dx
    Y = y + pieOffset[2] + pieR[2] * cos(phi) * dy
    lines(X, Y)

    # Add Center Point:
    points(x, y, col = "red", pch = 19, cex = 1.5)

    # Add Arrow:
    lines(c(x, x+pieOffset[1]), c(y, y+pieOffset[2]))

    # Add Color Wheel:
    psi = 2*pi*c(0, cumsum(abs(weights)/sum(abs(weights))))
    for (i in 1 : length(weights) ) {
        # Plotting Only Pie pieces with Weights > 5%
        if(psi[i+1]-psi[i] > 0.05 * 2*pi) {
            Psi = psi[i] + (0:100) * (psi[i+1]-psi[i])/100
            polyX = x + pieOffset[1] + pieR[1]*c(0, sin(Psi), 0) * dx
            polyY = y + pieOffset[2] + pieR[2]*c(0, cos(Psi), 0) * dy
            polygon(polyX, polyY, col = rainbow(nWeights)[i])
            # Adding the Asset Signs:
            text(x + pieOffset[1] + 0.75*pieR[1]* sin(Psi[51]) * dx,
                y + pieOffset[2] + 0.75*pieR[2]* cos(Psi[51]) * dy,
                col = "white", Sign[i])
         }
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.attributesWheel <-
    function(object, piePos = NULL, pieR = NULL, pieOffset = NULL, ...)
{
    # A function implemented by Oliver Greshake and Diethelm Wuertz

    # Description:
    #   Adds a pie plot of the weights

    # Arguments:
    
    # Details:
    #   The default settings are:
    #   piePos - Position of tangency Portfolio
    #   pieR - 10% of the Risk Range: diff(range(targetRisk(object)))/10

    # FUNCTION:

    # Extraction coordinates
    p = par()$usr/15
    dx = p[2]-p[1]
    dy = p[4]-p[3]

    # Pie Position:
    if(is.null(piePos)) {
        Data = getSeries(object)
        Spec = getSpec(object)
        Constraints = getConstraints(object)
        tg = getTargetReturn(tangencyPortfolio(Data, Spec, Constraints))
        ef = as.vector(getTargetReturn(object))
        piePos = which(diff(sign(ef-tg)) > 0)
    }

    # Pie Radius:
    if(is.null(pieR)) {
        pieR = c(1, 1)
    }

    # Pie Offset:
    if(is.null(pieOffset)) {
        pieOffset = c(2*dx, 0)
    }

    # Plot Circle - Get weighted Returns:
    weights = getWeights(object)
    dim = dim(weights)
    returns = getStatistics(object)$mu
    weightedReturns = NULL
    for(i in 1:dim[2]){
        nextWeightedReturns = weights[,i]*returns[i]
        weightedReturns = cbind(weightedReturns, nextWeightedReturns)
    }
    colnames(weightedReturns) = colnames(weights)
    weightedReturns = weightedReturns[piePos, ]
    nWeights = length(weightedReturns)
    Sign = rep("+", times =  nWeights)
    Sign[(1:nWeights)[weightedReturns < 0]] = "-"
    x = getTargetRisk(object)[piePos]
    y = getTargetReturn(object)[piePos]
    phi =  seq(0, 2*pi, length = 360)
    X = x + pieOffset[1] + pieR[1] * sin(phi) * dx
    Y = y + pieOffset[2] + pieR[2] * cos(phi) * dy
    lines(X, Y)

    # Add Center Point:
    points(x, y, col = "red", pch = 19, cex = 1.5)

    # Add Arrow:
    lines(c(x, x+pieOffset[1]), c(y, y+pieOffset[2]))

    # Add Color Wheel:
    psi = 2*pi*c(0, cumsum(abs(weightedReturns)/sum(abs(weightedReturns))))
    for (i in 1 : nWeights) {
        # Plotting Only Pie pieces with Weights > 5%
        if(psi[i+1]-psi[i] > 0.05 * 2*pi) {
            Psi = psi[i] + (0:100) * (psi[i+1]-psi[i])/100
            polyX = x + pieOffset[1] + pieR[1]*c(0, sin(Psi), 0) * dx
            polyY = y + pieOffset[2] + pieR[2]*c(0, cos(Psi), 0) * dy
            polygon(polyX, polyY, col = rainbow(nWeights)[i])
            # Adding the Asset Signs:
            text(x + pieOffset[1] + 0.75*pieR[1]* sin(Psi[51]) * dx,
                y + pieOffset[2] + 0.75*pieR[2]* cos(Psi[51]) * dy,
                col = "white", Sign[i])
         }
    }

    # Return Value:
    invisible()
}


#-------------------------------------------------------------------------------


.notStackedWeightsPlot <-
    function(object, col = NULL)
{
    # A function implemented by Oliver Greshake

    # Description:

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   col - a color palette, by default the rainbow palette

    # FUNCTION:

    # Settings:
    weights = getWeights(object)
    N = ncol(weights)
    targetRisk = getTargetRisk(object)[, 1]
    targetReturn = getTargetReturn(object)[, 1]
    nSigma = length(targetRisk)

    # Select Colors if not specified ...
    if (is.null(col)) col = rainbow(N)

    # Plot first asset ...
    plot(weights[, 1], col = col[1], type = "l", ylim = c(min(weights),
        max(weights)), xaxt = "n", xlab = "", ylab = "")

    # Add vertical Line at minimum risk:
    minIndex = which.min(targetRisk)
    minRisk = min(targetRisk)

    # Big Point at minimum risk for first asset ...
    points(x = minIndex, y = weights[minIndex, 1], col = col[1], pch = 19,
        xaxt = "n", yaxt = "n", cex = 2)

    # ... and all other assets
    for(i in 1:(N-1)){
        points(weights[, i+1], col = col[i+1], type = "l", xaxt = "n",
        yaxt = "n")
        points(x = minIndex, y = weights[minIndex, i+1], col = col[i+1],
            pch = 19, xaxt = "n", yaxt = "n", cex = 2)
    }
    grid()
    abline(h = 0, col = "grey", lty = 3)
    lines(x = c(minIndex, minIndex), y = c(0, 1), col = "black", lwd = 2)

    # Add Tailored Labels -  6 may be a good Number ...
    nLabels = 6
    M = c(0, ( 1: (nSigma %/% nLabels) ) ) * nLabels + 1
    text(minIndex, 1, "Min Risk", pos = 4)
    minRiskValue = as.character(signif(minRisk, 3))
    minReturnValue = as.character(signif(targetReturn[minIndex], 3))
    mtext(minRiskValue, side = 1, at = minIndex, cex = 0.7)
    mtext(minReturnValue, side = 3, line = 0.5, at = minIndex, cex = 0.7)

    # Take a reasonable number of significant digits to plot, e.g. 2 ...
    nPrecision = 3
    axis(1, at = M, labels = signif(targetRisk[M], nPrecision))
    axis(3, at = M, labels = signif(targetReturn[M], nPrecision))

    # Add Axis Labels and Title:
    mtext("Target Risk", side = 1, line = 2, cex = 0.7)
    mtext("Target Return", side = 3, line = 2, cex = 0.7)
    mtext("Weight", side = 2, line = 2, cex = 0.7)

    # Add Info:
    mtext(paste(getType(object), "|", getSolver(object)),
        side = 4, adj = 0, col = "grey", cex = 0.7)

    # Add Title:
    mtext("Weights", adj = 0, line = 2.5, font = 2, cex = 0.8)

    # Return Value:
    invisible()
}


#-------------------------------------------------------------------------------


.addlegend <-
    function(object, control = list())
{
    # A function implemented by Oliver Greshake

    # Description:
    #   Adds a perdefined legend to sliders

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   control - control list for colors and symbols

    # FUNCTION:

    # Settings:
    dim = getNAssets(object)
    namesSingleAsset = names(object@data$statistics$mu)
    # Check if polt is used for forntierSlider...
    if(control$sliderFlag == "frontier"){
        legendtext = c("Efficient Frontier", "Sharpe Ratio", "Minimum Variance",
            "Tangency Portfolio", "Market Portfolio", "Equal Weights",
            namesSingleAsset)
        color = c("black", control$sharpeRatio.col, control$minvariance.col,
            control$tangency.col, control$cml.col, control$equalWeights.col,
            control$singleAsset.col)
        sym = c(19, 19, control$minvariance.pch, control$tangency.pch,
            control$cml.pch, control$equalWeights.pch,
            rep(control$singleAsset.pch, times = dim))
    # ... else is the weightsSlider case
    } else {
            legendtext = c("Efficient Frontier", "Minimum Variance",
            "Tangency Portfolio", namesSingleAsset)
        color = c("black", control$minvariance.col,
            control$tangency.col, control$singleAsset.col)
        sym = c(19, control$minvariance.pch, control$tangency.pch,
            rep(control$singleAsset.pch, times = dim))
    }

    # Adding Legend:
    legend("topleft", legend = legendtext, col = color, pch = sym, cex = .8,
        bty = "n")

    # Return Value:
    invisible()
}


################################################################################

