
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
#  weightsPie                    Plots a pie of portfolio weights
#  weightedReturnsPie            Plots a pie of weighted means
#  covRiskBudgetsPie             Plots a pie of covariance risk budgets
#  tailRiskBudgetsPie            Plots a pie of copulae tail risk budgets
################################################################################


weightsPie <-
function(object, pos = NULL, labels = TRUE, col = NULL,
    box = TRUE, legend = TRUE, radius = 0.8, ...)
{
    # A function implemented by Diethelm Wuertz and Oliver Greshake

    # Description:
    #   Plots a Pie Chart of Weigths

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'.
    #   pos - a numeric value, determining the position on the efficient
    #       frontier plotting the pie, by default NULL, i.e. expecting
    #       an object having only one set of weights like the tangency
    #       portfolio.
    #   box - a logical value, determining whether a frame (box) should
    #       be plotted around the pie, by default TRUE.
    #   col - a color palette, by default the rainbow palette.
    #   legend - a logical value, determining whether a legend with
    #       the names of the assets should be plotted, by default TRUE.

    # Example:
    #   weightsPie(tangencyPortfolio(dutchPortfolioData(), portfolioSpec()))
    #   title(main = "Tangency Portfolio Weights")

    # FUNCTION:

    # Default Settings:
    Title <- "Weights"
    if (is.null(col)) col <- seqPalette(getNAssets(object), "Blues")
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX = 0.9 else CEX = 0.7

    # Get Weights:
    if (is.null(pos)) {
        Weights <- getWeights(object@portfolio)
    } else {
        Weights <- getWeights(object@portfolio)[pos, ]   
    }  
    X <- Weights
    
    # Check for Negative Pie Segments:
    nX <- getNAssets(object)
    Sign <- rep("+", nX)
    Sign[(1:nX)[X < 0]] <- "-"
    absX <- abs(X)
    Index <- (1:nX)[X > 0]

    # Take care of labels, they are also used by the function pie():
    if (!is.logical(labels)) {
        Names <- pieLabels <- labels
        labels <- FALSE
    } else  {
        Names <- pieLabels <- object@data@data$names
    }

    # Pie Chart:
    col <- col[Index]
    legendAssets <- Names[Index]
    Labels <- paste(Names, Sign)
    Labels = Labels[X > 0]
    Y <- X[X > 0]

    # Plot:
    if (labels) {
        pie(Y, labels = Labels, col = col, radius = radius, cex = CEX)
    } else {
        pie(Y, labels = pieLabels, col = col, radius = radius, ...)
    }

    # Add Title:
    if (labels) 
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)

    # Add Info:
    if (labels) {
        mtext(paste(getType(object), "|", getSolver(object)),
            side = 4, adj = 0, col = "grey", cex = 0.7)
    }

    # Add Legend:
    if (legend) {
        legend("topleft", legend = legendAssets, bty = "n", cex = CEX,
            fill = col)
        legendY <- as.character(round(100*Y, digits = 1))
        legendY <- paste(Sign[Index], legendY, sep = "")
        legendY <- paste(legendY, "%")
        legend("topright", legend = legendY, bty = "n", cex = CEX,
            fill = col)
    }

    # Add Box:
    if (box) box()

    # Return Value:
    invisible(Y)
}


# ------------------------------------------------------------------------------


weightedReturnsPie <-
function(object, pos = NULL, labels = TRUE, col = NULL,
    box = TRUE, legend = TRUE, radius = 0.8, ...)
{
    # A function implemented by Diethelm Wuertz and Oliver Greshake

    # Description:
    #   Adds a pie plot of the weights

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'.
    #   pos - a numeric value, determining the position on the efficient
    #       frontier plotting the pie, by default NULL, i.e. expecting
    #       an object having only one set of weights like the tangency
    #       portfolio.
    #   box - a logical value, determining whether a frame (box) should
    #       be plotted around the pie, by default TRUE.
    #   col - a color palette, by default the rainbow palette.
    #   legend - a logical value, determining whether a legend with
    #       the names of the assets should be plotted, by default TRUE.

    # Example:
    #   attributesPie(tangencyPortfolio(dutchPortfolioData(), portfolioSpec()))
    #   title(main = "Tangency Portfolio Weights")

    # FUNCTION:

    # Default Settings:
    Title <- "Weighted Returns"
    if (is.null(col)) col <- seqPalette(getNAssets(object), "Blues")
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX = 0.9 else CEX = 0.7


    # Get Weights:
    if (is.null(pos)) {
        Weights <- getWeights(object@portfolio)
    } else {
        Weights <- getWeights(object@portfolio)[pos, ]   
    }  
    Returns = getStatistics(object)$mu
    X <- Weights * Returns

    # Check for Negative Pie Segments:
    nX <- getNAssets(object)
    Sign <- rep("+", nX)
    Sign[(1:nX)[X < 0]] <- "-"
    absX <- abs(X)
    Index <- (1:nX)[X > 0]

    # Take care of labels, they are also used by the function pie():
    if (!is.logical(labels)) {
        Names <- pieLabels <- labels
        labels <- FALSE
    } else  {
        Names <- pieLabels <- object@data@data$names
    }

    # Pie Chart:
    col <- col[Index]
    legendAssets <- Names[Index]
    Labels <- paste(Names, Sign)
    Labels <- Labels[X > 0]
    Y <- X[X > 0]

    # Plot:
    if (labels) {
        pie(Y, labels = Labels, col = col, radius = radius, cex = CEX)
    } else {
        pie(Y, labels = pieLabels, col = col, radius = radius, ...)
    }

    # Add Title:
    if (labels) 
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)

    # Add Info:
    if (labels) {
        mtext(paste(getType(object), "|", getSolver(object)),
            side = 4, adj = 0, col = "grey", cex = 0.7)
    }

    # Add Legend:
    if (legend) {
        legend("topleft", legend = legendAssets, bty = "n", cex = CEX, 
        	fill = col)
        legendY = as.character(round(100*Y, digits = 1))
        legendY = paste(Sign[Index], legendY, sep = "")
        legendY = paste(legendY, "%")
        legend("topright", legend = legendY, bty = "n", cex = CEX,
            fill = col)
    }

    # Add Box:
    if (box) box()

    # Return Value:
    invisible(Y)
}


# ------------------------------------------------------------------------------


covRiskBudgetsPie <-
function(object, pos = NULL, labels = TRUE, col = NULL,
    box = TRUE, legend = TRUE, radius = 0.8, ...)
{
    # A function implemented by Diethelm Wuertz and Oliver Greshake

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'.
    #   pos - a numeric value, determining the position on the efficient
    #       frontier plotting the pie, by default NULL, i.e. expecting
    #       an object having only one set of weights like the tangency
    #       portfolio.
    #   box - a logical value, determining whether a frame (box) should
    #       be plotted around the pie, by default TRUE.
    #   col - a color palette, by default the rainbow palette.
    #   legend - a logical value, determining whether a legend with
    #       the names of the assets should be plotted, by default TRUE.

    # Description:
    #   Plots a Pie Chart of Risk Budgets

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   col - a color palette, by default the rainbow palette

    # Example:
    #   riskBudgetsPie(tangencyPortfolio(dutchPortfolioData(), portfolioSpec()))
    #   title(main = "Tangency Portfolio Weights")

    # FUNCTION:

    # Default Settings:
    Title = "Covariance Risk Budgets"
    if (is.null(col)) col <- seqPalette(getNAssets(object), "Blues")
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX <- 0.9 else CEX <- 0.7

    # Get Cov Risk Budgets:
    if (is.null(pos)) {
        X <- getCovRiskBudgets(object@portfolio)
    } else {
        X <- getCovRiskBudgets(object@portfolio)[pos, ]   
    }  

    # Check for Negative Pie Segments:
    nX <- getNAssets(object)
    Sign <- rep("+", nX)
    Sign[(1:nX)[X < 0]] <- "-"
    absX <- abs(X)
    Index <- (1:nX)[X > 0]

    # Take care of labels, they are also used by the function pie():
    if (!is.logical(labels)) {
        Names <- pieLabels <- labels
        labels <- FALSE
    } else  {
        Names <- pieLabels <- object@data@data$names
    }

    # Legend Labels:
    col <- col[Index]
    legendAssets <- Names[Index]
    Labels <- paste(Names, Sign)
    Labels <- Labels[X > 0]
    Y <- X[X > 0]

    # Plot:
    if (labels) {
        pie(Y, labels = Labels, col = col, radius = radius, cex = CEX)
    } else {
        pie(Y, labels = pieLabels, col = col, radius = radius, ...)
    }

    # Add Title:
    if (labels) 
         mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)

    # Add Info:
    if (labels) {
        mtext(paste(getType(object), "|", getSolver(object)),
            side = 4, adj = 0, col = "grey", cex = 0.7)
    }

    # Add Legend:
    if (legend) {
        legend("topleft", legend = legendAssets, bty = "n", cex = CEX,
            fill = col)
        legendY <- as.character(round(100*Y, digits = 1))
        legendY <- paste(Sign[Index], legendY, sep = "")
        legendY <- paste(legendY, "%")
        legend("topright", legend = legendY, bty = "n", cex = CEX,
            fill = col)
    }

    # Add Box:
    if (box) box()

    # Return Value:
    invisible(Y)
}


# ------------------------------------------------------------------------------


tailRiskBudgetsPie <-
function(object, pos = NULL, labels = TRUE, col = NULL,
    box = TRUE, legend = TRUE, radius = 0.8, ...)
{
    ### todo: take care of @portfolio slot ...
 
    # A function implemented by Diethelm Wuertz and Oliver Greshake

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'.
    #   pos - a numeric value, determining the position on the efficient
    #       frontier plotting the pie, by default NULL, i.e. expecting
    #       an object having only one set of weights like the tangency
    #       portfolio.
    #   box - a logical value, determining whether a frame (box) should
    #       be plotted around the pie, by default TRUE.
    #   col - a color palette, by default the rainbow palette.
    #   legend - a logical value, determining whether a legend with
    #       the names of the assets should be plotted, by default TRUE.

    # Description:
    #   Plots a Pie Chart of Tail Risk Budgets

    # Arguments:
    #   object - an object of class 'fPORTFOLIO'
    #   col - a color palette, by default the rainbow palette

    # Example:
    #   riskBudgetsPie(tangencyPortfolio(dutchPortfolioData(), portfolioSpec()))
    #   title(main = "Tangency Portfolio Weights")

    # FUNCTION:

    # Default Settings:
    Title <- "Tail Risk Budgets"
    if (is.null(col)) col <- seqPalette(getNAssets(object), "Blues")
    if (sum(c(par()$mfrow, par()$mfcol)) == 4) CEX <- 0.9 else CEX <- 0.7

    # Extracting weights position, if specified
    if(!is.null(pos)){
        object = object
        object@portfolio$weights = getWeights(object@portfolio)[pos, ]
    }

    # Check:
    stop("Not yet implemented")
    tailRiskMatrix = getTailRisk(object)
    X <- getCovRiskBudgets(object)

    # Check for Negative Pie Segments:
    nX <- getNAssets(object)
    Sign <- rep("+", nX)
    Sign[(1:nX)[X < 0]] <- "-"
    absX <- abs(X)
    Index <- (1:nX)[X > 0]

    # Take care of labels, they are also used by the function pie():
    if (!is.logical(labels)) {
        Names <- pieLabels <- labels
        labels <- FALSE
    } else  {
        Names <- pieLabels <- object@data@data$names    
    }

    # Legend Labels:
    col <- col[Index]
    legendAssets <- Names[Index]
    Labels <- paste(Names, Sign)
    Labels <- Labels[X > 0]
    Y <- X[X > 0]

    # Plot:
    if (labels) {
        pie(Y, labels = Labels, col = col, radius = radius, cex = CEX)
    } else {
        pie(Y, labels = pieLabels, col = col, radius = radius, ...)
    }

    # Add Title:
    if (labels) 
        mtext(Title, adj = 0, line = 2.5, font = 2, cex = CEX+0.1)

    # Add Info:
    if (labels) {
        mtext(paste(getType(object), "|", getSolver(object)),
            side = 4, adj = 0, col = "grey", cex = 0.7)
    }

    # Add Legend:
    if (legend) {
        legend("topleft", legend = legendAssets, bty = "n", cex = CEX,
            fill = col)
        legendY = as.character(round(100*Y, digits = 1))
        legendY = paste(Sign[Index], legendY, sep = "")
        legendY = paste(legendY, "%")
        legend("topright", legend = legendY, bty = "n", cex = CEX,
            fill = col)
    }

    # Add Box:
    if (box) box()

    # Return Value:
    invisible(Y)
}


################################################################################

