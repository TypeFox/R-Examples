
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
# FUNCTION:              DESCRIPTION:   
#  .emaIndicator           Exponential moving average indicator
#  .macdIndicator          MACD Indicator
#  .drawdownsIndicator     Maximum drawdowns Indicator
# FUNCTION:              DESCRIPTION:   
#  .rebalancingStats       Rebalancing statistics
############################################################################### 


.emaIndicator <- 
    function(series, lambda)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Exponential moving average indicator    
    
    # FUNCTION:
    
    # EMA:
    x <- rep(mean(series[1:10,]), times=nrow(series))
    for (i in 2:nrow(series))
        x[i] <- (1-lambda)*series[i] + lambda*x[i-1]
    x <- as.timeSeries(data=x, charvec=time(series), units=colnames(series))
    
    # Return Value:
    x
}


# -----------------------------------------------------------------------------


.macdIndicator <- 
    function(index, spar=0.5, lambda=c(0.80, 0.85, 0.90), 
        trace = TRUE, doplot=TRUE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   MACD price/index indicator
    
    # FUNCTION:
    
    # Series:
    rets <- returns(index)
    Index <- log(index)[-1, ]
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, 
        main = "MACD Analytics", trace=TRUE, doplot=FALSE)
    ablines <- tps$ablines
    
    # MACD Analytics:
    ema1 <- .emaIndicator(Index, lambda=lambda[1])
    ema2 <- .emaIndicator(Index, lambda=lambda[2])
    macd <- ema1 - ema2
    signal <- .emaIndicator(macd, lambda[3])
    histogram <- macd - signal
    
    # Indicator:
    indicator <- sign(histogram)
    rebalancing <- .rebalancingStats(index, indicator, trace=trace)
      
    # Plot Turning Points:
    if(doplot) {
        turnsAnalytics(index=index, spar=spar, 
            main="MACD Index Indicator", 
            trace=FALSE, doplot=doplot)
        tradePositions <- as.vector(indicator)
        tradeForecasts <- c(0, tradePositions[-length(tradePositions)])
        outSample <- Index[1] + log(cumulated(rets*tradeForecasts))
        Ups <- Index[as.vector(indicator) == 1, ]
        if(nrow(Ups) > 0) points(Ups, pch=19, cex=0.33, col="green")   
        Downs <- Index[as.vector(indicator) == 0, ]
        if(nrow(Downs) > 0) points(Downs, pch=19, cex=0.33, col="blue")   
        lines(outSample, col="magenta")
        box(col="white")
        box(bty="l")
    }
    
    # Plot Indicator:
    if(doplot) {
        plot(macd, col="green", ylab=paste("MACD", colnames(index)))
        abline(v=ablines, lty=3, lwd=2, col="grey")
        lines(histogram, type="h", col="black")
        lines(macd, col="red")
        mtext(paste("lambda: ", lambda[1], lambda[2], lambda[2], sep=" "), 
            adj=0, side=4, cex=0.7, col="darkgrey")    
        box(col="white")
        box(bty="l") 
    }
    
    # Return Value:
    invisible(list(index=index, macd=macd, histogram=histogram, 
        rebalancing=rebalancing))
}


# -----------------------------------------------------------------------------


.drawdownsIndicator <- 
    function(index, spar=0.5, lambda=c(0.80, 0.85, 0.10),  
        trace=TRUE, doplot=TRUE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Drawdown analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    
    # FUNCTION:
    
    # Series:
    rets <- returns(index)
    Index <- log(index)[-1, ]
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, 
        main = "MACD Drawdown Analytics", trace=TRUE, doplot=FALSE)
    ablines <- tps$ablines
    
    # Returns and Drawdowns:
    dd <- drawdowns(rets)
       
    # Long/Short Drawdowns EMA:
    mdd1 <- .emaIndicator(dd, lambda[1])
    mdd2 <- .emaIndicator(dd, lambda[2])
        
    # MACD/Signal/Histogram Line:
    macd <- mdd1 - mdd2
    signal <- .emaIndicator(macd, lambda[3])
    histogram <- macd - signal
        
    # Indicator:
    indicator <- rets
    series(indicator) <- 1-sign(as.integer(macd < 0 & histogram < 0 )) 
    rebalancing <- .rebalancingStats(index, indicator, trace=trace)
    
    # Plot Turning Points:
    if(doplot) {
        tps <- turnsAnalytics(index=index, spar=spar, 
            main="MACD Drawdown Indicator", trace=FALSE, doplot=doplot)
        tradePositions <- as.vector(indicator)
        tradeForecasts <- c(0, tradePositions[-length(tradePositions)])
        outSample <- Index[1] + log(cumulated(rets*tradeForecasts))
        Ups <- Index[as.vector(indicator) == 1,]
        if(nrow(Ups) > 0) points(Ups, pch=19, cex=0.33, col="green")   
        Downs <- Index[as.vector(indicator) == 0, ]
        if(nrow(Downs) > 0) points(Downs, pch=19, cex=0.33, col="blue")   
        lines(outSample, col="magenta")
        box(col="white")
        box(bty="l")
    }
       
    # Plot Indicator:
    if(doplot) {
        plot(mdd1, ylim=c(min(dd), max(macd)), ylab=colnames(index))
        positions <- tps$positions
        ablines <- tps$ablines
        abline(v=ablines, lty=1, lwd=2, col="lightgrey")    
        Time <- time(indicator)
        Time <- Time[!as.logical(indicator)] 
        abline(v=Time, lty=3, lwd=2, col="steelblue")         
        lines(mdd1, col="black")
        lines(mdd2, col="red")
        lines(max(abs(macd))*histogram/max(abs(histogram)), type="h", col="orange")      
        lines(max(abs(macd))*histogram/max(abs(histogram)), type="l", col="orange")      
        abline(h=0, col="grey") 
        mtext(paste("lambda: ", lambda[1], lambda[2], lambda[2], sep=" "), 
            adj=0, side=4, cex=0.7, col="darkgrey")            
        box(col="white")
        box(bty="l")
    }
    
    # Return Value:
    invisible(list(indicator=indicator, index=index, returns=rets, 
        drawdowns=dd, macd=macd, signal=signal, histogram=histogram, 
        rebalancing=rebalancing))
}


###############################################################################    


.rebalancingStats <- 
    function(index, indicator, trace=TRUE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Simple rebalancing statistics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    
    # FUNCTION:
    
    # Returns:
    rets <- returns(index)
    
    # Rebalancing:
    tradePositions <- as.vector(indicator)
    tradeForecasts <- c(0, tradePositions[-length(tradePositions)])
    rebalancing <- c(
        max=sum(abs(rets)),
        insample=sum(rets*tradePositions),
        forecasts=sum(rets*tradeForecasts),
        rets=sum(rets))      
        
    # Trace:
    if (trace) {
        cat("Rebalancing:\n")
        print(rebalancing)
    }
    
    # Return Value:
    invisible(rebalancing)
}


############################################################################### 

