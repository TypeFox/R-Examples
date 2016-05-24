
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
# FUNCTION:                     DESCRIPTION:                
#  weightsLinePlot               Plots staggered weights
#  weightedReturnsLinePlot       Plots staggered weighted returns
#  covRiskBudgetsLinePlot        Plots covariance risk budgets
################################################################################


weightsLinePlot <- 
function(object, labels = TRUE, col = NULL, title = TRUE, 
    box = TRUE, legend = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots a bar chart of weights
    
    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   labels - should the graph be automatically labeled?
    #   col - a color palette, by default the rainbow palette
    #   title - should the graph get default title and labels?
    #   legend - should a legend be added to the plot?
    
    # FUNCTION:
    
    # Use default color if not specified ...
    Title = "Weights"
    if (is.null(col)) col = seqPalette(getNAssets(object)+1, "Blues")[-1]
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX = 0.9 else CEX = 0.7
    
    # Compute Weights:
    X = weights = getWeights(object)
    
    # Define Plot Range:
    ymax = max(colMaxs(weights))
    ymin = min(colMins(weights))
    range = ymax - ymin
    ymax = ymax + 0.005 * range
    ymin = ymin - 0.005 * range
    dim = dim(weights)
    range = dim[1]
    xmin = 0
    xmax = range + 0.2 * range
    
    # Create Bar Plots:
    if (labels) {
        if(legend){
            ts.plot(X, 
                gpars = list(col = col, ann = FALSE, xaxt = "n"), 
                xlim = c(xmin, xmax), ylim = c(ymin, ymax))
            legendtext = names(getStatistics(object)$mu)
            if(is.null(legendtext)){
                for(i in 1:dim[2]){legendtext[i] = paste("Asset", i, sep = " ")}
            }
            legend("topright", legend = legendtext, bty = "n", cex = CEX,
                fill = col)
        } else {
            ts.plot(weights, gpars = list(col = col, ann = FALSE, xaxt = "n"))
        }
    } else {
        ts.plot(X, ...)
    }
    
    # Add Tailored Labels -  6 may be a good Number ...
    targetRisk = getTargetRisk(object)[, 1]
    targetReturn = getTargetReturn(object)[, 1]
    nSigma = length(targetRisk)
    nLabels = 6
    M = c(0, ( 1:(nSigma %/% nLabels) ) ) *nLabels + 1
    nSignifDigits = 3
    axis(3, at = M, labels = signif(targetRisk[M], nSignifDigits))
    axis(1, at = M, labels = signif(targetReturn[M], nSignifDigits))
    
    # Add Axis Labels and Title:
    if (title) {
        mtext("Target Risk", side = 3, line = 2, adj = 1, cex = CEX)
        mtext("Target Return", side = 1, line = 2, adj = 1, cex = CEX)
        mtext("Weight", side = 2, line = 2, adj = 1, cex = CEX)
    }
      
    # Add Weights 0 and 1 Reference Lines
    # lines(x = c(0, nSigma), c(1, 1), col = "grey", lty = 3) 
    # lines(x = c(0, nSigma), c(0, 0), col = "grey", lty = 3)   
    
    # Add vertical Line at minimum risk:
    minIndex = which.min(targetRisk)
    minRisk = signif(min(targetRisk), 3)
    abline(v = minIndex, col = "black", lty = 1, lwd = 2)
    
    # Add Info:
    if (title) {
        mtext(paste(getType(object), "|", getSolver(object), "|", "minRisk =", 
            minRisk), side = 4, adj = 0, col = "grey", cex = 0.7)
    }
        
    # Add Title:
    if (title) {
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)
    }
    
    # Complete to draw box ...
    if (box) box()
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


weightedReturnsLinePlot <- 
function(object, labels = TRUE, col = NULL, title = TRUE, 
    box = TRUE, legend = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots a bar chart of weights
    
    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   labels - should the graph be automatically labeled?
    #   col - a color palette, by default the rainbow palette
    #   title - should the graph get default title and labels?
    #   legend - should a legend be added to the plot?
    
    # FUNCTION:
    
    # Use default color if not specified ...
    Title = "Weighted Returns"
    if (is.null(col)) col = seqPalette(getNAssets(object)+1, "Blues")[-1]
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX = 0.9 else CEX = 0.7
    
    # Compute Weighted Returns:
    weights = getWeights(object)
    dim = dim(weights)
    returns = getStatistics(object)$mu
    weightedReturns = NULL
    for(i in 1:dim[2]){
        nextWeightedReturns = weights[,i]*returns[i]
        weightedReturns = cbind(weightedReturns, nextWeightedReturns)
    }
    colnames(weightedReturns) = colnames(weights)
    X = weightedReturns
    
    # Define Plot Range:
    ymax = max(colMaxs(X))
    ymin = min(colMins(X))
    range = ymax - ymin
    ymax = ymax + 0.005 * range
    ymin = ymin - 0.005 * range
    dim = dim(weights)
    range = dim[1]
    xmin = 0
    xmax = range + 0.2 * range
    
    # Create Bar Plots:
    if (labels) {
        if(legend){
            ts.plot(X, 
                gpars = list(col = col, ann = FALSE, xaxt = "n"), 
                xlim = c(xmin, xmax), ylim = c(ymin, ymax))
            legendtext = names(getStatistics(object)$mu)
            if(is.null(legendtext)){
                for(i in 1:dim[2]){legendtext[i] = paste("Asset", i, sep = " ")}
            }
            legend("topright", legend = legendtext, bty = "n", cex = CEX,
                fill = col)
        } else {
            ts.plot(weights, gpars = list(col = col, ann = FALSE, xaxt = "n"))
        }
    } else {
        ts.plot(X, ...)
    }
    
    # Add Tailored Labels -  6 may be a good Number ...
    targetRisk = getTargetRisk(object)[, 1]
    targetReturn = getTargetReturn(object)[, 1]
    nSigma = length(targetRisk)
    nLabels = 6
    M = c(0, ( 1:(nSigma %/% nLabels) ) ) *nLabels + 1
    nSignifDigits = 3
    axis(3, at = M, labels = signif(targetRisk[M], nSignifDigits))
    axis(1, at = M, labels = signif(targetReturn[M], nSignifDigits))
    
    # Add Axis Labels and Title:
    if (title) {
        mtext("Target Risk", side = 3, line = 2, adj = 1, cex = CEX)
        mtext("Target Return", side = 1, line = 2, adj = 1, cex = CEX)
        mtext("Weight", side = 2, line = 2, adj = 1, cex = CEX)
    }
      
    # Add Weights 0 and 1 Reference Lines
    # lines(x = c(0, nSigma), c(1, 1), col = "grey", lty = 3) 
    # lines(x = c(0, nSigma), c(0, 0), col = "grey", lty = 3)   
    
    # Add vertical Line at minimum risk:
    minIndex = which.min(targetRisk)
    minRisk = signif(min(targetRisk), 3)
    abline(v = minIndex, col = "black", lty = 1, lwd = 2)
    
    # Add Info:
    if (title) {
        mtext(paste(getType(object), "|", getSolver(object), "|", "minRisk =", 
            minRisk), side = 4, adj = 0, col = "grey", cex = 0.7)
    }
        
    # Add Title:
    if (title) {
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)
    }
    
    # Complete to draw box ...
    if (box) box()
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


covRiskBudgetsLinePlot <- 
function(object, labels = TRUE, col = NULL, title = TRUE, 
    box = TRUE, legend = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots a bar chart of weights
    
    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   labels - should the graph be automatically labeled?
    #   col - a color palette, by default the rainbow palette
    #   title - should the graph get default title and labels?
    #   legend - should a legend be added to the plot?
    
    # FUNCTION:
    
    # Use default color if not specified ...
    Title = "Covariance Risk Budgets"
    if (is.null(col)) col = seqPalette(getNAssets(object)+1, "Blues")[-1]
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX = 0.9 else CEX = 0.7
    
    # Compute Covariance Risk Budgets:
    X = getCovRiskBudgets(object)
 
    # Define Plot Range:
    ymax = max(colMaxs(X))
    ymin = min(colMins(X))
    range = ymax - ymin
    ymax = ymax + 0.005 * range
    ymin = ymin - 0.005 * range
    dim = dim(X)
    range = dim[1]
    xmin = 0
    xmax = range + 0.2 * range
    
    # Create Bar Plots:
    if (labels) {
        if(legend){
            ts.plot(X, 
                gpars = list(col = col, ann = FALSE, xaxt = "n"), 
                xlim = c(xmin, xmax), ylim = c(ymin, ymax))
            legendtext = names(getStatistics(object)$mu)
            if(is.null(legendtext)){
                for(i in 1:dim[2]){legendtext[i] = paste("Asset", i, sep = " ")}
            }
            legend("topright", legend = legendtext, bty = "n", cex = CEX,
                fill = col)
        } else {
            ts.plot(weights, gpars = list(col = col, ann = FALSE, xaxt = "n"))
        }
    } else {
        ts.plot(X, ...)
    }
    
    # Add Tailored Labels -  6 may be a good Number ...
    targetRisk = getTargetRisk(object) 
    targetReturn = getTargetReturn(object) 
    nSigma = length(targetRisk)
    nLabels = 6
    M = c(0, ( 1:(nSigma %/% nLabels) ) ) *nLabels + 1
    nSignifDigits = 3
    axis(3, at = M, labels = signif(targetRisk[M], nSignifDigits))
    axis(1, at = M, labels = signif(targetReturn[M], nSignifDigits))
    
    # Add Axis Labels and Title:
    if (title) {
        mtext("Target Risk", side = 3, line = 2, adj = 1, cex = CEX)
        mtext("Target Return", side = 1, line = 2, adj = 1, cex = CEX)
        mtext("Weight", side = 2, line = 2, adj = 1, cex = CEX)
    }
      
    # Add Weights 0 and 1 Reference Lines
    # lines(x = c(0, nSigma), c(1, 1), col = "grey", lty = 3) 
    # lines(x = c(0, nSigma), c(0, 0), col = "grey", lty = 3)   
    
    # Add vertical Line at minimum risk:
    minIndex = which.min(targetRisk[, 1])
    minRisk = signif(min(targetRisk[, 1]), 3)
    abline(v = minIndex, col = "black", lty = 1, lwd = 2)
    
    # Add Info:
    if (title) {
        mtext(paste(getType(object), "|", getSolver(object), "|", "minRisk =", 
            minRisk), side = 4, adj = 0, col = "grey", cex = 0.7)
    }
        
    # Add Title:
    if (title) {
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)
    }
    
    # Complete to draw box ...
    if (box) box()
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.tailRiskBudgetsLinePlot <- 
function(object, labels = TRUE, col = NULL, title = TRUE, 
    box = TRUE, legend = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots a bar chart of weights
    
    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   labels - should the graph be automatically labeled?
    #   col - a color palette, by default the rainbow palette
    #   title - should the graph get default title and labels?
    #   legend - should a legend be added to the plot?
    
    # FUNCTION:
    
    # Use default color if not specified ...
    Title = "Covariance Risk Budgets"
    if (is.null(col)) col = seqPalette(getNAssets(object)+1, "Blues")[-1]
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX = 0.9 else CEX = 0.7
    
    # Compute Covariance Risk Budgets:
    stop("Not yet implemented")
    tailRiskMatrix = getTailRisk(object) 
    X = getTailRiskBudgets(object)
 
    # Define Plot Range:
    ymax = max(colMaxs(X))
    ymin = min(colMins(X))
    range = ymax - ymin
    ymax = ymax + 0.005 * range
    ymin = ymin - 0.005 * range
    dim = dim(weights)
    range = dim[1]
    xmin = 0
    xmax = range + 0.2 * range
    
    # Create Bar Plots:
    if (labels) {
        if(legend){
            ts.plot(X, 
                gpars = list(col = col, ann = FALSE, xaxt = "n"), 
                xlim = c(xmin, xmax), ylim = c(ymin, ymax))
            legendtext = names(getStatistics(object)$mu)
            if(is.null(legendtext)){
                for(i in 1:dim[2]){legendtext[i] = paste("Asset", i, sep = " ")}
            }
            legend("topright", legend = legendtext, bty = "n", cex = CEX,
                fill = col)
        } else {
            ts.plot(weights, gpars = list(col = col, ann = FALSE, xaxt = "n"))
        }
    } else {
        ts.plot(X, ...)
    }
    
    # Add Tailored Labels -  6 may be a good Number ...
    targetRisk = getTargetRisk(object) 
    targetReturn = getTargetReturn(object) 
    nSigma = length(targetRisk)
    nLabels = 6
    M = c(0, ( 1:(nSigma %/% nLabels) ) ) *nLabels + 1
    nSignifDigits = 3
    axis(3, at = M, labels = signif(targetRisk[M], nSignifDigits))
    axis(1, at = M, labels = signif(targetReturn[M], nSignifDigits))
    
    # Add Axis Labels and Title:
    if (title) {
        mtext("Target Risk", side = 3, line = 2, adj = 1, cex = CEX)
        mtext("Target Return", side = 1, line = 2, adj = 1, cex = CEX)
        mtext("Weight", side = 2, line = 2, adj = 1, cex = CEX)
    }
      
    # Add Weights 0 and 1 Reference Lines
    # lines(x = c(0, nSigma), c(1, 1), col = "grey", lty = 3) 
    # lines(x = c(0, nSigma), c(0, 0), col = "grey", lty = 3)   
    
    # Add vertical Line at minimum risk:
    minIndex = which.min(targetRisk[, 1])
    minRisk = signif(min(targetRisk[, 1]), 3)
    abline(v = minIndex, col = "black", lty = 1, lwd = 2)
    
    # Add Info:
    if (title) {
        mtext(paste(getType(object), "|", getSolver(object), "|", "minRisk =", 
            minRisk), side = 4, adj = 0, col = "grey", cex = 0.7)
    }
        
    # Add Title:
    if (title) {
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)
    }
    
    # Complete to draw box ...
    if (box) box()
    
    # Return Value:
    invisible()
}


################################################################################

