
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


###############################################################################
# FUNCTION:                DESCRIPTION:
#  markowitzHull            Hull for a long-only Markowitz portfolio
#  feasibleGrid             Square grid on top of the feasible set
#  bestDiversification      Diversified portfolios on top of the feasible set
#  riskSurface              Risk values on top of the feasible set
#  surfacePlot              Risk Suface Plot for a Markowitz portfolio
# FUNCTION:                DESCRIPTION:
#  .scaledColors            Quantile color scaling 
################################################################################


markowitzHull <-
    function (data, nFrontierPoints=50)
{
    # Description:
    #    Returns the Hull for a long-only Markowitz portfolio
    
    # Arguments:
    #    data - an object of class 'timeSeries'
   
    # Example:
    #    hull <- markowitzHull(100*LPP2005.RET[, 1:6], nFrontierPoints=11) 
    #    plot(hull[[1]], type="n"); polygon(hull[[1]], col="grey")
    
    # FUNCTION:
    
    # Check:
    stopifnot(is.timeSeries(data))   
    
    # Compute Frontier and Minimum Variance Locus:
    Spec <- portfolioSpec()
    setNFrontierPoints(Spec) <- nFrontierPoints
    frontier <- portfolioFrontier(data, spec=Spec)
    Risks <- risks <- frontierPoints(frontier)[, 1]
    Returns <- frontierPoints(frontier)[, 2]
    
    # Compute Maximum Variance Locus - Pairwise Assets Approach:
    N <- ncol(data)
    for (i in 1:(N - 1)) for (j in (i + 1):N) {
        Data <- data[, c(i, j)]
        ans <- portfolioFrontier(Data, spec=Spec)
        coord <- frontierPoints(ans)
        nextFrontier <- approx(coord[, 2], coord[, 1], xout = Returns)$y
        naIndex <- which(is.na(nextFrontier))
        nextFrontier[naIndex] <- Risks[naIndex]
        risks <- rbind(risks, nextFrontier)
    }
    
    # Hull:
    targetReturn <- Returns
    minTargetRisk <- Risks
    maxTargetRisk <- colMaxs(risks)
    hull <- cbind(
        targetReturn = Returns, 
        minTargetRisk = Risks,
        maxTargetRisk = colMaxs(risks))
    
    # Polygon:
    polygon <- cbind(
        c(minTargetRisk, rev(maxTargetRisk)[-1]), 
        c(targetReturn, rev(targetReturn)[-1]) )
    rownames(polygon) <- 1:nrow(polygon)
    colnames(polygon) <- c("targetRisk", "targetReturn")
    
    # Return Value:
    ans <- polygon 
    attr(ans, "data") <- data
    attr(ans, "hull") <- hull
    attr(ans, "frontier") <- frontier
    invisible(ans)
}


# -----------------------------------------------------------------------------


feasibleGrid <-
    function(hull, trace=FALSE)
{
    # Description:
    #    Returns best diversified portfolios on top of the feasible Set
    
    # Arguments:
    #    hull - an object as returned from the function markowitzHull
    #    trace - a logical, should the function be traced ?
    
    # Example:
    #    hull <- markowitzHull(100*LPP2005.RET[, 1:6], nFrontierPoints=21) 
    #    grid <- feasibleGrid(hull, TRUE) 
    
    # FUNCTION:
    
    # Data:
    polygon <- hull
    data <- attr(hull, "data")
    hull <- attr(hull, "hull")
    
    # Trace Hull:
    if (trace) {
        plot(polygon)
        box(col="white")
        polygon(polygon, col="grey")
        grid()
    }
    
    # Settings:
    minRisks <- as.vector(hull[, 2])
    maxRisks <- as.vector(hull[, 3])
    minRisk <- min(minRisks)
    maxRisk <- max(maxRisks)
    targetRisks <- seq(minRisk, maxRisk, length = length(minRisks))
    targetReturns <- as.vector(hull[, 1])
    N <- length(targetReturns)
    
    # Get Weights on Grid:
    Grid <- matrix(NA, ncol=N, nrow=N)
    offset <- diff(range(targetRisks[1:2]))/2
    for (i in 1:N) {
        targetReturn <- targetReturns[i]
        for (j in 1:N) {
            targetRisk <- targetRisks[j] + offset
            if (targetRisk >= minRisks[i] && targetRisk <= maxRisks[i]) {
                Grid[j, i] <- 1
                if (trace) points(targetRisk, targetReturn, pch=19)
            }
       }
    }
    
    # Return Value:
    ans <- list(x=targetRisks, y=targetReturns, z=Grid)
    attr(ans, "data") <- data
    attr(ans, "polygon") <- polygon
    attr(ans, "hull") <- hull
    class(ans) <- c("feasibleGrid", "list")
    invisible(ans)
}


# -----------------------------------------------------------------------------


bestDiversification <-
    function(grid, FUN="var", trace=FALSE)
{
    # Description:
    #    Returns best diversified portfolios on top of the feasible Set
    
    # Arguments:
    #    data - an object of class 'timeSeries'
    #    grid - an object of class 'feasibleGrid' 
    #        as returned by the function feasibleGid()
    #    FUN - the divesification function, a function with
    #        with the weights as its first argument
    #    trace - a logical, should the function be traced ?
    
    # Example:
    #    data <- 100*LPP2005.RET[, 1:6]
    #    hull <- makowitzHull(data, nFrontierPoints=21) 
    #    grid <- feasibleGrid(hull, trace=TRUE) 
    #    diversification <- bestDiversification(grid, FUN=var)
    
    # FUNCTION:
    
    # Data:
    data <- attr(grid, "data")
    polygon <- attr(grid, "polygon")
       
    # Settings:
    targetRisks <- grid$x
    targetReturns <- grid$y
    Grid <- grid$z
    N <- length(targetRisks)
    objectiveFun <- match.fun(FUN)
    nAssets <- ncol(data)
    MEAN <- colMeans(data)
    COV <- cov(data)
    
    # Trace:
    if(trace) {
       image(grid, col="lightgrey")
       box(col="white")
       grid()
    }
    
    # Get Weights on Grid:
    Weights <- Coord <- NULL
    Objective <- NA * Grid
    Start <- rep(1/nAssets, times = nAssets)
    for (i in 1:N) {
        targetReturn <- targetReturns[i]
        for (j in 1:N) {
            targetRisk <- targetRisks[j]
            if (!is.na(Grid[j,i])) {
                ans <- donlp2NLP(
                    start = Start,
                    objective <- objectiveFun,
                    par.lower = rep(0, times = nAssets), 
                    par.upper = rep(1, times = nAssets),
                    eqA = rbind(rep(1, times = nAssets), MEAN),
                    eqA.bound = c(1, targetReturn),
                    eqFun = list(function(x) sqrt(t(x) %*% COV %*% x)),
                    eqFun.bound = targetRisk)
                Weights <- rbind(Weights, ans$solution)
                Objective[j,i] <- objectiveFun(ans$solution)
                Coord <- rbind(Coord, c(j,i))
                if(trace) {
                    points(targetRisk, targetReturn, pch=19, cex=0.7)
                }
            }
        }
    }
   
    # Return Value:
    ans <- list(x=targetRisks, y=targetReturns, z=Objective)
    attr(ans, "data") <- data
    attr(ans, "polygon") <- polygon
    attr(ans, "weights") <- cbind(Coord, Weights)
    class(ans) <- c("bestDiversification", "list")
    invisible(ans)
}


# -----------------------------------------------------------------------------


riskSurface <- 
    function(diversification, FUN=NULL, ...)
{
    # Description:
    #    Returns a risk values on top of the feasible set
    
    # Arguments:
    #    diversification - an object of class class 'bestDiversification'
    #       as returned by the function bestDiversification()
    #    FUN - risk surface function having arguments
    #       FUN(data, weights, ...)
    #    ... - optional arguments passed to FUN
    
    # Example:
    #    data <- 100*LPP2005.RET[, 1:6]
    #    hull <- markowitzHull(data, nFrontierPoints=21) 
    #    grid <- feasibleGrid(hull, TRUE) 
    #    diversification <- bestDiversification(grid)
    #    surface <- riskSurface(diversification)
    
    # FUNCTION:
    
    # Data and Weighs:
    data <- attr(diversification, "data")
    weights <- attr(diversification, "weights")
    polygon <- attr(diversification, "polygon")
    
    # Risk Function:
    if (is.null(FUN)) FUN <- function(data, weights, ...) var(weights)
    fun <- match.fun(FUN)
    
    # Grid:
    Coord <- attr(diversification, "weights")[, 1:2]
    Weights <- attr(diversification, "weights")[, -(1:2)]
    N <- nrow(Coord)
    x <- diversification$x
    y <- diversification$y
    z <- diversification$z
    
    # Risk Surface:
    Value <- NA * z
    for (k in 1:N) {
        Value[Coord[k, 1], Coord[k, 2]] <- fun(data, Weights[k, ], ...) 
    }
    
    # Return Value:
    ans <- list(x=x, y=y, z=Value)
    attr(ans, "data") <- data
    attr(ans, "weights") <- weights
    attr(ans, "polygon") <- polygon
    class(ans) <- c("riskSurface", "list")
    ans
}


###############################################################################


surfacePlot <- 
    function(surface, type=c("image", "filled.contour"),
    nlevels=11, palette=topo.colors, addContour=TRUE, addGrid=TRUE,
    addHull=TRUE, addAssets=TRUE, ...)
{
    # Description:
    
    # Arguments:
    #    surface - an object of class 'riskSurface' as 
    #        returned by the function riskSurface()
    #    type - a character string denoting the plot type,
    #        by default "image, alternatively "filledContour"
    #    nlevels - integer, the number of countour levels
    #    palette - color palette function
    #    addCountour - a logical flag, should contour lines be added ?
    #    addCountour - a logical flag, should contour lines be added ?
    #    addGrid - a logical flag, should grid lines be added ?
    #    addAssets - a logical flag, should assets points be added ?
    #    ... - optional arguments passed to the function title()
    
    # Example:
    #    data <- 100*LPP2005.RET[, 1:6]
    #    hull <- markowitzHull(data) 
    #    grid <- feasibleGrid(hull, trace=TRUE) 
    #    diversification <- bestDiversification(grid, trace=TRUE)
    #    surface <- riskSurface(diversification)
    #    surfacePlot(surface)
    #    surfacePlot(surface, type="f"); 
    #         title("Weights Diversification", xlab="Risk", ylab="Return")
    
    # FUNCTION:
    
    # Surface Points;
    x <- surface$x
    y <- surface$y
    z <- surface$z

    # Quantile Levels:
    colors <- .scaledColors(surface, palette=palette, nlevels=nlevels)
    levels <- colors$levels
    palette <- colors$palette

    # Contour overlayed Image Ranges:  
    yOffset <- 0.025*diff(range(y))
    yLim <- c(min(y)-yOffset, max(y)+yOffset)
    xOffset <- 0.1*diff(range(x))
    xLim <- c(min(x)-xOffset/4, max(x)+xOffset)

    # Select Type:
    type <- match.arg(type)
    if (type == "image") {
        image(x, y, z, xlim=xLim, ylim=yLim, xlab="", ylab="", col=palette)
        box(col="white")
    }
    else if (type == "filled.contour") {
        image(x, y, z, xlim=xLim, ylim=yLim, xlab="", ylab="", col="white")
        
        # DW
        # .Internal(filledcontour()) no longer works on 3.0.
        # .Internal(filledcontour(
        #     as.double(x), as.double(y), z, 
        #     as.double(levels), col = palette))
        
        # Use instead:
        graphics::.filled.contour(
          x = as.double(x), 
          y = as.double(y), 
          z = z,
          levels = as.double(levels), 
          col = palette)
        
        box(col="white")
    }
    
    # Add Contour Lines:
    if(addContour) contour(x, y, z, add=TRUE, levels=signif(levels, 3))
    
    # Add Hull:
    if(addHull) {
      hull <- attr(surface, "polygon")
      lines(hull, lwd=2, col="darkgreen")
    }
    
    # Add Grid:
    if(addGrid) grid()
    
    # Add Optional Lables:
    title(...)

    # Add Legend:
    cs <- cumsum(levels)
    css <- ( cs - min(cs) ) / diff(range(cs))
    css <- 0.95 * css + 0.025
    cy <- min(y) + css * diff(range(y))
    cx <- rep(xLim[2]-0.1 * xOffset, length(cy))
    lines(cx, cy, lwd=3)
    for (i in 1:(nlevels-1)) lines(c(cx[i], cx[i+1]), c(cy[i], cy[i+1]), lwd=3, col=palette[i])
    for (i in 1:nlevels) points(cx[i], cy[i], pch=16, cex=1.1, col="black")
    textOffset <- c(-0.0005, 0.0005, 0.0008, 0.0008, rep(0, 7))
    text(cx, cy+textOffset, as.character(signif(levels, 2)), pos=2, cex=0.8)
      
    # Add Assets:
    if (addAssets) {
        frontier <- portfolioFrontier(data)
        pointCex <- 2.5
        textCex <- 0.5
        xy <- minvariancePoints(frontier, auto=FALSE, 
           pch=19, cex=pointCex, col = "red")
        text(xy[, 1], xy[, 2], "MVP", font=2, col="white", cex=textCex)
        xy <- tangencyPoints(frontier, auto=FALSE, 
           pch=19, cex=pointCex, col="orange")
        text(xy[, 1], xy[, 2], "TGP", font=2, col="white", cex=textCex)
        xy <- equalWeightsPoints(frontier, auto=FALSE, 
           pch=19, cex=pointCex, col="brown")
        text(xy[, 1], xy[, 2], "EWP", font=2, col="white", cex=textCex)
        xy <- singleAssetPoints(frontier, auto=FALSE, 
           pch=19, cex=pointCex, col="black", lwd=2)
        text(xy[, 1], xy[, 2], rownames(xy), font=2, col="white", cex=textCex)
    }
       
    # Return Value:
    invisible(list(surface=surface, levels=levels))
}


# ------------------------------------------------------------------------------


.scaledColors <- 
    function(surface, palette=topo.colors, nlevels=11)
{
    # Description:
    #   scales a color palette
    
    # Arguments:
    #    surface - a list with x,y positions and z values
    #    palette - color palette function
    #    bin - quantile bin width of contour levels
    
    # FUNCTION:
    
    # Extract Surface Risk Values:
    Z <- as.vector(surface$z)
  
    # Scale by Equidistant Quantiles:
    levels <- quantile(Z, probs=seq(from=0, to=1, length=nlevels), na.rm=TRUE)
    
    # Compose Color Palette:
    palette <- palette(nlevels-1)
  
    # Return Value:
    list(palette=palette, levels=levels)
}


###############################################################################


