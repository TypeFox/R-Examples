
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


###############################################################################
# FUNCTION:            DESCRIPTION:
#  ternaryMap           Displays a risk map for ternary portfolios     
#  ternaryFrontier      Plots the efficient frontier of a ternary portfolio
# FUNCTION:            DESCRIPTION:
#  riskMap              normalVaR risk map function called from ternaryMap()
#  maxddMap             max Drawdown risk map function called from ternaryMap()
# FUNCTION:            DESCRIPTION:
#  ternaryWeights       Creates a set of ternary weights
#  ternaryCoord         Computes x, y coordinates from weights
#  ternaryPoints        Adds points to a ternary map plot
###############################################################################


ternaryMap <- 
    function(data, FUN=NULL, ..., 
    locator=FALSE, N=41, palette=topo.colors, nlevels=11)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Displays a risk map for ternary portfolios
    
    # Arguments:
    #    data - a ternary 'timeSeries' object of financial returns
    #    FUN - the map function
    #    ... - optional arguments passed to function FUN
    #    locator - a logical flag to activate the locator
    #    N - number of bins
    #    palette - color palette
    #    nlevels - number of contour levels
    
    # Examples:
    #    ternaryMap(data=SWX.RET[, 1:3])
    #    ternaryMap(data=SWX.RET[, 1:3], FUN=.riskTest, locator=TRUE)
    #    ternaryMap(data=SWX.RET[, 1:3], FUN=.maxddTest)
    
    # FUNCTION:
    
    # N=41; palette=topo.colors; nlevels=11; FUN <- NULL; ... <- NULL
    
    # Surface Function:
    if (is.null(FUN)) FUN <- function(data, weights, ...) var(weights)
    fun <- match.fun(FUN)
    
    # Grid Points:
    s <- sqrt(3)/2
    x <- seq(0, 1, length=N)
    y <- s * x 
    
    # Surface Values:
    G <- matrix(NA, N, N)
    xy <- W <- NULL
    for (j in 1:N) {
        w3 <- y[j] / s
        for (i in 1:N) {
            w2 <- x[i] - w3/2
            w1 <- 1 - w2 - w3
            if (w1 >= -1/N && w2 >= -1/N) {
            # if (w1 >= 0 && w2 >= 0) {
               w <- c(w1, w2, w3)
               xy <- rbind(xy, c(x[i], y[j]))
               W <- rbind(W, w)
               G[i, j] <- ans <- fun(data, w, ...)
            }
        }
    }
    surface <- list(x=x, y=y/s, z=G)   
    x <- surface$x
    y <- surface$y
    z <- surface$z

    # Color Settings:
    colors <- .scaledColors(surface, palette = palette, nlevels = nlevels)
    levels <- colors$levels
    palette <- colors$palette
    
    # Image Ranges:
    yOffset <- 0.025 * diff(range(y))
    yLim <- c(min(y) - yOffset, max(y) + yOffset)
    xOffset <- 0.1 * diff(range(x))
    xLim <- c(min(x) - xOffset/4, max(x) + xOffset)
    
    # Filled Contour Plot:
    image(x, y, z, xlim = xLim, ylim = yLim, xlab = "", ylab = "", 
          col = "white")
    grid()
    
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
    
    contour(x, y, z, add = TRUE, levels = signif(levels, 3))
    d <- 3/N
    polygon(c(1, 1+d, 0.5+d, 0.5, 1), c(0, 0, s, s, 0)/s, col="white", 
            border="white")
    polygon(c(0, 0.5, 0.5-d,  -d, 0), c(0, s, s, 0, 0)/s, col="white", 
            border="white")
    box(col = "white")
    
    # Please do not Remove:
    mtext("Rmetrics", 4, col="grey", adj=0, cex=0.7)

    # Grid Lines:
    for(k in 1:10) 
      lines(c(k*0.1, k*0.05), c(0, k*0.1), col="grey", lty=3)
    for(k in 1:10) 
      lines(c(1-k*0.1, 1-k*0.05), c(0, k*0.1), col="grey", lty=3)
    for(k in 1:9) 
      lines(c(k*0.05, 1-k*0.05), c(k*0.1, k*0.1), col="grey", lty=3)

    # Add Legend:
    cs <- cumsum(levels)
    css <- (cs - min(cs))/diff(range(cs))
    css <- 0.95 * css + 0.025
    cy <- min(y) + css * diff(range(y))
    cx <- rep(xLim[2] - 0.1 * xOffset, length(cy))
    lines(cx, cy, lwd = 3)
    for (i in 1:(nlevels - 1)) lines(c(cx[i], cx[i + 1]), c(cy[i], 
        cy[i + 1]), lwd = 3, col = palette[i])
    for (i in 1:nlevels) points(cx[i], cy[i], pch = 16, cex = 1.1, 
        col = "black")
    textOffset <- c(-5e-04, 5e-04, 8e-04, 8e-04, rep(0, 7))
    text(cx, cy + textOffset, as.character(signif(levels, 2)), 
        pos = 2, cex = 0.8)

    # Decoration:
    pointCex <- 2.5
    textCex <- 0.5
    
    # EfficientFrontier:
    frontier <- portfolioFrontier(data)
    weights <- getWeights(frontier@portfolio)
    xy <- cbind(x=weights[, 2]+weights[, 3]/2, y=weights[, 3])
    lines(xy, col="brown", lwd=3)
    
    # Global Minimum Variance Portfolio:
    weights<- getWeights(minvariancePortfolio(data))
    xy <- c(x=weights[2]+weights[3]/2, y=weights[3])
    points(xy[1], xy[2], pch=19, cex=pointCex, col="red")
    text(xy[1], xy[2], "MVP", font=2, col="white",  cex=textCex)
    
    # Tangency Portfolio:
    weights <- getWeights(tangencyPortfolio(data))
    xy <- c(x=weights[2]+weights[3]/2, y=weights[3])
    points(xy[1], xy[2], pch=19, cex=pointCex, col="orange")
    text(xy[1], xy[2], "TGP", font=2, col="white", cex=textCex)
    
    # Equal Weights:
    xy <- rbind(c(1/3+1/6, 1/3))
    points(xy, pch=19, cex=pointCex, col="brown")
    text(xy, "EWP", font=2, col="white", cex = textCex)
     
    # Individual Assets:   
    xy = rbind(c(0,0), c(1,0), c(1/2, 1))
    points(xy, pch=19, cex=pointCex, col="black")
    text(xy, colnames(data), font = 2, col = "white", cex = textCex)

    # Locator:
    if (locator) {
        for (i in 1:512) {
            ans <- locator(n = 1, type = "p", pch=10) 
            w3 <- ans$y
            w2 <- ans$x - w3/2
            w <- round(c(w1=1-w2-w3, w2, w3), 2)
            names(w) <- colnames(data)
            total <- sum(w)
            names(total) <- "Total"
            z <- signif(fun(data, w, ...), 3)
            ans <- data.frame(rbind(c(w, total, z)))
            rownames(ans) <- "Composition"
            print(ans)
        }
    }
    
    # Return Value:
    invisible()
}


# -----------------------------------------------------------------------------


ternaryFrontier <-
    function(data, locator=FALSE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Plots the efficient frontier of a ternary map
    
    # Arguments:
    #    data - a ternary 'timeSeries' object of financial returns
    #    locator - a logical flag to activate the locator
    
    # Example:
    #    ternaryFrontier(SWX.RET[, 1:3], locator=TRUE)
    
    # FUNCTION:
    
    # Long Only Markowitz Portfolio:
    polygon <- markowitzHull(data)
    object <- attr(polygon, "frontier")
      
    # Plot Range:
    offset <- 0.1
    xlim <- c(0, max(sqrt(diag(getCov(object)))))
    Xlim <- c(xlim[1] - diff(xlim) * offset, xlim[2] + diff(xlim) * offset)
    ylim <- range(getMean(object))
    Ylim <- c(ylim[1] - diff(ylim) * offset, ylim[2] + diff(ylim) * offset)

    # Get Points and Add Frontier:
    
    frontierPlot(object, auto=FALSE, xlim=Xlim, ylim=Ylim, pch=19, 
                 labels=FALSE)
    polygon(polygon, col="grey", border="grey")
    points(frontierPoints(object, frontier="upper"), pch=19)
    box(col="white")
    grid()
    
    # Please do not Remove:
    mtext("Rmetrics", 4, col="grey", adj=0, cex=0.7)

    # Zero  Axis Lines:
    abline(h = 0, col = "grey")
    abline(v = 0, col = "grey")
           
    # Decoration Settings:
    pointCex <- 2.5
    textCex <- 0.5

    # Add Global Minimum Variance Portfolio:
    xy <- minvariancePoints(
      object, auto=FALSE, pch=19, col="red", cex=pointCex)
    text(xy[1], xy[2], "MVP", font=2, col="white",  cex=textCex)
      
    # Add Tangency Portfolio for zero risk free Rate:
    tangencyLines(object, auto=FALSE, col="blue")  
    xy <- tangencyPoints(
      object, auto=FALSE, pch=19, col="blue", cex=pointCex) 
    text(xy[1], xy[2], "TGP", font=2, col="white",  cex=textCex)
    
    # Add Equal Weights Portfolio:   
    xy <- equalWeightsPoints(
      object, auto=FALSE, pch=19, col="brown", cex=pointCex)
    text(xy[1], xy[2], "EWP", font=2, col="white", cex=textCex)
    
    # Add Two Assets Portfolios: 
    xy <- singleAssetPoints(
      object, auto=FALSE, pch=19, col="black", cex=pointCex)
    text(xy[,1], xy[,2], colnames(data), font=2, col="white", cex=textCex)

    # Add Sharpe Ratio:
    sharpeRatioLines(object, auto=FALSE, col="orange", lwd=2, pch=19)
    
    # Locator:
    if (locator) {
        for (i in 1:512) {
            ans <- locator(n = 1, type = "p", pch=10) 
            Risk <- ans$x
            Return <- ans$y
            SharpeRatio <- ans$y/ans$x
            ans <- data.frame(rbind(c(Risk, Return, SharpeRatio)))
            colnames(ans) <- c("Risk", "Return", "SharpeRatio")
            rownames(ans) <- "Portfolio"
            print(signif(ans, 3))
        }
    }
     
    # Retirn Value:
    invisible(object)
}


# #############################################################################


riskMap <- 
    function(data, weights) 
{  
    # Description:
    #    normalVaR risk map function called from ternaryMap()
    
    # FUNCTION:
    
    # Map:
    tS <- pfolioReturn(data, weights=as.vector(weights))
    ans <- normalVaR(tS) 
    names(ans) <- "normalVaR"
    
    # Return Value:
    ans
}


# -----------------------------------------------------------------------------


maxddMap <- 
    function(data, weights) 
{  
    # Description:
    #    Max Drawdown map function called from ternaryMap()
    
    # FUNCTION:
    
    # Map:
    tS <- pfolioReturn(data, weights=as.vector(weights))
    ans <- colMins(tS) 
    names(ans) <- "maxdd"
    
    # Return Value:
    ans
}
    

# #############################################################################
# Utility Functions:


ternaryWeights <- 
    function(n=21)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns a set of ternary weights
    
    # Arguments:
    #    n - number of bins
    
    # Example:
    #    ternaryWeights()
    
    # FUNCTION:
    
    # Settings:
    eps <- sqrt(.Machine$double.eps)
    
    # Creates a set of ternary weights
    W <- seq(0, 1, length=n)
    W1 <- rep(W, times = length(W))
    W2 <- rep(W, each = length(W))
    W3 <- 1 - W1 - W2
    W3[abs(W3) < eps] <- 0
    W <- cbind(W1, W2, W3)
    weights <- W[W1 + W2 <= 1, ]
    
    # Return Value:
    weights
}


# -----------------------------------------------------------------------------


ternaryCoord <- 
    function(weights)  
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns x,y Coordinates of weights triangle
    
    # Example:
    #    ternaryCoord(ternaryWeights())
    
    # FUNCTION:
    
    # Computexs x, y coordinates from weights
    x <- 1 - weights[,1] - weights[, 2]/2
    y <- sqrt(3) * weights[, 2] /2
    
    # Return Value:
    cbind(x=x, y=y)
}


# -----------------------------------------------------------------------------

  
ternaryPoints <- 
    function(weights, ...)  
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Adds points to a ternary map plot
    
    # Example:
    #    ternaryPoints(ternaryWeights())
    
    # FUNCTION:
    
    # Transpose in the case of a single weight
    if(is.null(dim(weights))) weights <- t(weights)
    
    # Adds points to a ternary map
    coord <- ternaryCoord(weights)
    points(coord, ...)
    
    # Return Value:
    invisible()
}


###############################################################################


