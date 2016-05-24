### R code from vignette source 'timeSeriesPlot.Rnw'

###################################################
### code chunk number 1: environment
###################################################
Sys.setlocale("LC_ALL", "C")


###################################################
### code chunk number 2: library
###################################################
require(timeSeries)
require(xts)
require(PerformanceAnalytics)
require(fTrading)
tS1 <- 100 * cumulated(LPP2005REC[, 1])     # SBI (univariate)
tS2 <- 100 * cumulated(LPP2005REC[, 1:2])   # SBI & SPI (bivariate)
tS3 <- 100 * cumulated(LPP2005REC[, 1:3])   # SBI, SPI, SWIIT (Swiss Market)
tS6 <- 100 * cumulated(LPP2005REC[, 1:6])   # Swiss and Foreign Market Indexes


###################################################
### code chunk number 3: univariateSingle
###################################################
par(mfrow=c(1, 1))
plot(tS1) 


###################################################
### code chunk number 4: univariateSinglePlot
###################################################
par(mfrow=c(1, 1))
plot(tS1) 


###################################################
### code chunk number 5: univariateSingle2
###################################################
require(PerformanceAnalytics)
par(mfrow=c(3, 1))
xts::plot.xts(as.xts(tS1)) 
PerformanceAnalytics::chart.TimeSeries(as.xts(tS1)) 
plot(tS1) 


###################################################
### code chunk number 6: univariateSingle2Plot
###################################################
require(PerformanceAnalytics)
par(mfrow=c(3, 1))
xts::plot.xts(as.xts(tS1)) 
PerformanceAnalytics::chart.TimeSeries(as.xts(tS1)) 
plot(tS1) 


###################################################
### code chunk number 7: multivariateSingle
###################################################
par(mfrow=c(1, 1))             
plot(tS3, plot.type="s")   


###################################################
### code chunk number 8: multivariateSinglePlot
###################################################
par(mfrow=c(1, 1))             
plot(tS3, plot.type="s")   


###################################################
### code chunk number 9: multivariateSingle2
###################################################
par(mfrow=c(2, 1))       
require(PerformanceAnalytics)
PerformanceAnalytics::chart.TimeSeries(as.xts(tS3))   
plot(tS3, plot.type="s")   


###################################################
### code chunk number 10: multivariateSingle2Plot
###################################################
par(mfrow=c(2, 1))       
require(PerformanceAnalytics)
PerformanceAnalytics::chart.TimeSeries(as.xts(tS3))   
plot(tS3, plot.type="s")   


###################################################
### code chunk number 11: oneColMultiple
###################################################
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m")  


###################################################
### code chunk number 12: oneColMultiplePlot
###################################################
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m")  


###################################################
### code chunk number 13: twoColMultiple
###################################################
par(mfrow=c(1, 1))
plot(tS6, plot.type="m")  


###################################################
### code chunk number 14: twoColMultiplePlot
###################################################
par(mfrow=c(1, 1))
plot(tS6, plot.type="m")  


###################################################
### code chunk number 15: gapMultiple
###################################################
par(mfrow=c(1, 1))
plot(tS3, plot.type="m", mar=c(gap=0.3, 5.1, gap=0.3, 2.1)) 


###################################################
### code chunk number 16: gapMultiplePlot
###################################################
par(mfrow=c(1, 1))
plot(tS3, plot.type="m", mar=c(gap=0.3, 5.1, gap=0.3, 2.1)) 


###################################################
### code chunk number 17: combineSingle
###################################################
par(mfrow=c(2, 1))
par(mar = c(bottom=1.5, 5.1, top=4, 2.1))
plot(tS2[, 1])
par(mar = c(bottom=4, 5.1, top=1.5, 2.1))
plot(tS2[, 2])


###################################################
### code chunk number 18: combineSinglePlot
###################################################
par(mfrow=c(2, 1))
par(mar = c(bottom=1.5, 5.1, top=4, 2.1))
plot(tS2[, 1])
par(mar = c(bottom=4, 5.1, top=1.5, 2.1))
plot(tS2[, 2])


###################################################
### code chunk number 19: layoutSingle
###################################################
nf <- layout(mat=matrix(c(1, 1, 2, 3), byrow = TRUE, nrow=2))
par(mar = c(bottom=2, 5.1, top=3, 2.1))
plot(tS3[, 1])
par(mar = c(bottom=3, 5.1, top=2, 1.1))
plot(tS3[, 2])
par(mar = c(bottom=3, 4.1, top=2, 2.1))
plot(tS3[, 3])


###################################################
### code chunk number 20: layoutSinglePlot
###################################################
nf <- layout(mat=matrix(c(1, 1, 2, 3), byrow = TRUE, nrow=2))
par(mar = c(bottom=2, 5.1, top=3, 2.1))
plot(tS3[, 1])
par(mar = c(bottom=3, 5.1, top=2, 1.1))
plot(tS3[, 2])
par(mar = c(bottom=3, 4.1, top=2, 2.1))
plot(tS3[, 3])


###################################################
### code chunk number 21: layout2Single
###################################################
nf <- layout(mat=matrix(c(1, 1, 2, 3), byrow=TRUE, nrow=2), heights=c(2.5,1))
par(mar = c(bottom=2, 5.1, top=3, 2.1))
plot(tS3[, 1])
par(mar = c(bottom=3, 5.1, top=1.5, 1.1))
plot(tS3[, 2])
par(mar = c(bottom=3, 4.1, top=1.5, 2.1))
plot(tS3[, 3])


###################################################
### code chunk number 22: layout2SinglePlot
###################################################
nf <- layout(mat=matrix(c(1, 1, 2, 3), byrow=TRUE, nrow=2), heights=c(2.5,1))
par(mar = c(bottom=2, 5.1, top=3, 2.1))
plot(tS3[, 1])
par(mar = c(bottom=3, 5.1, top=1.5, 1.1))
plot(tS3[, 2])
par(mar = c(bottom=3, 4.1, top=1.5, 2.1))
plot(tS3[, 3])


###################################################
### code chunk number 23: scatter
###################################################
par(mfrow=c(1,1))
plot(tS2[, 1], tS2[, 2])


###################################################
### code chunk number 24: scatterPlot
###################################################
par(mfrow=c(1,1))
plot(tS2[, 1], tS2[, 2])


###################################################
### code chunk number 25: pretty
###################################################
par(mfcol = c(2, 1))
plot(tS1, at = "pretty")
plot(tS1, at = "chic")  


###################################################
### code chunk number 26: prettyPlot
###################################################
par(mfcol = c(2, 1))
plot(tS1, at = "pretty")
plot(tS1, at = "chic")  


###################################################
### code chunk number 27: chicUnivariateSingle
###################################################
par(mfcol=c(2, 1))
plot(tS3, plot.type="s", at="pretty")
plot(tS3, plot.type="s", at="chic")


###################################################
### code chunk number 28: chicUnivariateSinglePlot
###################################################
par(mfcol=c(2, 1))
plot(tS3, plot.type="s", at="pretty")
plot(tS3, plot.type="s", at="chic")


###################################################
### code chunk number 29: minorTicks
###################################################
par(mfrow=c(3, 1))                
plot(tS1, minor.ticks="day", at="pretty")
plot(tS1, minor.ticks="week", at="pretty")
plot(tS1, minor.ticks="month", at="pretty")


###################################################
### code chunk number 30: minorTicksPlot
###################################################
par(mfrow=c(3, 1))                
plot(tS1, minor.ticks="day", at="pretty")
plot(tS1, minor.ticks="week", at="pretty")
plot(tS1, minor.ticks="month", at="pretty")


###################################################
### code chunk number 31: chicOneColMultiple
###################################################
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m", at="pretty")


###################################################
### code chunk number 32: chicOneColMultiplePlot
###################################################
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m", at="pretty")


###################################################
### code chunk number 33: chicTwoColMultiple
###################################################
par(mfrow=c(1, 1)) 
plot(tS6, plot.type="m", at="chic")  


###################################################
### code chunk number 34: chicTwoColMultiplePlot
###################################################
par(mfrow=c(1, 1)) 
plot(tS6, plot.type="m", at="chic")  


###################################################
### code chunk number 35: tailoredAxis
###################################################
par(mfrow=c(2, 1))
at <- paste0("200", c("6-01", "6-04", "6-07", "6-10", "7-01", "7-04"), "-01")
plot(tS3, plot.type="s", format="%B\n%Y", at=at)
plot(tS3, plot.type="s", format="%b/%y", at=at)


###################################################
### code chunk number 36: tailoredAxisPlot
###################################################
par(mfrow=c(2, 1))
at <- paste0("200", c("6-01", "6-04", "6-07", "6-10", "7-01", "7-04"), "-01")
plot(tS3, plot.type="s", format="%B\n%Y", at=at)
plot(tS3, plot.type="s", format="%b/%y", at=at)


###################################################
### code chunk number 37: annSingle
###################################################
par(mfrow=c(2, 2))
plot(tS1, ann=FALSE)                
plot(tS3, plot.type="s", ann=FALSE, at="pretty")                                   
plot(tS6, plot.type="s", ann=FALSE, at="pretty")


###################################################
### code chunk number 38: annSinglePlot
###################################################
par(mfrow=c(2, 2))
plot(tS1, ann=FALSE)                
plot(tS3, plot.type="s", ann=FALSE, at="pretty")                                   
plot(tS6, plot.type="s", ann=FALSE, at="pretty")


###################################################
### code chunk number 39: titleSingle
###################################################
par(mfrow=c(2, 2))
plot(tS1); title(main = "Index") 
plot(tS3, plot.type="s"); title(main = "Index") 
plot(tS3, plot.type="s"); title(main = "Index", xlab = "Date")                              
plot(tS6, plot.type="s"); title(main = "Index", xlab = "Date")  


###################################################
### code chunk number 40: titleSinglePlot
###################################################
par(mfrow=c(2, 2))
plot(tS1); title(main = "Index") 
plot(tS3, plot.type="s"); title(main = "Index") 
plot(tS3, plot.type="s"); title(main = "Index", xlab = "Date")                              
plot(tS6, plot.type="s"); title(main = "Index", xlab = "Date")  


###################################################
### code chunk number 41: axisFontSize
###################################################
par(mfrow=c(3, 1))
plot(tS3, at="chic", plot.type="s", cex.axis=0.75)                            
plot(tS3, at="chic", plot.type="s", cex.axis=1.00)                                   
plot(tS3, at="chic", plot.type="s", cex.axis=1.25)


###################################################
### code chunk number 42: axisFontSizePlot
###################################################
par(mfrow=c(3, 1))
plot(tS3, at="chic", plot.type="s", cex.axis=0.75)                            
plot(tS3, at="chic", plot.type="s", cex.axis=1.00)                                   
plot(tS3, at="chic", plot.type="s", cex.axis=1.25)


###################################################
### code chunk number 43: flipAxisOne
###################################################
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m", yax.flip = TRUE)                                    


###################################################
### code chunk number 44: flipAxisOnePlot
###################################################
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m", yax.flip = TRUE)                                    


###################################################
### code chunk number 45: typeMultiple
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", type=c("l", "p", "h"), at="pretty")


###################################################
### code chunk number 46: typeMultiplePlot
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", type=c("l", "p", "h"), at="pretty")


###################################################
### code chunk number 47: colorNamesMultiple
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", col=c("blue", "orange", "darkgreen"))


###################################################
### code chunk number 48: colorNamesMultiplePlot
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", col=c("blue", "orange", "darkgreen"))


###################################################
### code chunk number 49: palettesMultiple
###################################################
par(mfrow=c(1, 1))                                                
plot(tS6, plot.type="s", col=heat.colors(n=6, alpha = 1), 
  at="chic", format = "%B\n%Y")


###################################################
### code chunk number 50: palettesMultiplePlot
###################################################
par(mfrow=c(1, 1))                                                
plot(tS6, plot.type="s", col=heat.colors(n=6, alpha = 1), 
  at="chic", format = "%B\n%Y")


###################################################
### code chunk number 51: ltyMultiple
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", col=1, lty=1:3, at="chic")


###################################################
### code chunk number 52: ltyMultiplePlot
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", col=1, lty=1:3, at="chic")


###################################################
### code chunk number 53: lwdMultiple
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", col=1, lwd=3:1, at="chic")


###################################################
### code chunk number 54: lwdMultiplePlot
###################################################
par(mfrow=c(1, 1))                                                
plot(tS3, plot.type="m", col=1, lwd=3:1, at="chic")


###################################################
### code chunk number 55: symbolsSizeMultiple
###################################################
par(mfrow=c(1, 1))    
plot(tS3, plot.type="s", type="p", 
  col=1:3, pch=21:23, cex.pch=c(0.2, 0.2, 0.2), at="pretty")


###################################################
### code chunk number 56: symbolsSizeMultiplePlot
###################################################
par(mfrow=c(1, 1))    
plot(tS3, plot.type="s", type="p", 
  col=1:3, pch=21:23, cex.pch=c(0.2, 0.2, 0.2), at="pretty")


###################################################
### code chunk number 57: gridSingle
###################################################
par(mfrow=c(1, 1))    
plot(tS3, plot.type="s", grid=FALSE)


###################################################
### code chunk number 58: gridSinglePlot
###################################################
par(mfrow=c(1, 1))    
plot(tS3, plot.type="s", grid=FALSE)


###################################################
### code chunk number 59: noBoxSingle
###################################################
par(mfrow=c(1, 1))    
plot(tS3, plot.type="s", frame.plot=FALSE, grid=FALSE)
box()
box(bty = "7", col = "white") # boxL
grid(NA, NULL, col = "darkgrey") # hgrid


###################################################
### code chunk number 60: gridSinglePlot
###################################################
par(mfrow=c(1, 1))    
plot(tS3, plot.type="s", grid=FALSE)


###################################################
### code chunk number 61: horizMultiple
###################################################
par(mfrow=c(1, 1))
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, col=col)
  abline(h=0, col = "brown", lwd=2)}
plot(returns(tS3), plot.type="m", col = .colorwheelPalette(3),
  panel=lines2, at="pretty")


###################################################
### code chunk number 62: horizMultiplePlot
###################################################
par(mfrow=c(1, 1))
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, col=col)
  abline(h=0, col = "brown", lwd=2)}
plot(returns(tS3), plot.type="m", col = .colorwheelPalette(3),
  panel=lines2, at="pretty")


###################################################
### code chunk number 63: rugMultiple
###################################################
par(mfrow=c(1, 1))
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, type="h", col=col)
  rug(Y, side=4, col="steelblue") }
plot(returns(tS6), plot.type="m", col = .colorwheelPalette(6), 
  panel=lines2, at="pretty")


###################################################
### code chunk number 64: rugMultiplePlot
###################################################
par(mfrow=c(1, 1))
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, type="h", col=col)
  rug(Y, side=4, col="steelblue") }
plot(returns(tS6), plot.type="m", col = .colorwheelPalette(6), 
  panel=lines2, at="pretty")


###################################################
### code chunk number 65: emaMultiple
###################################################
par(mfrow=c(1, 1))
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, type="l", col=col)
  lines(x=X, y=emaTA(Y), col="black") }
plot(tS3, plot.type="m", col = .colorwheelPalette(3), panel=lines2, 
  grid=TRUE, at="pretty")


###################################################
### code chunk number 66: emaMultiplePlot
###################################################
par(mfrow=c(1, 1))
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, type="l", col=col)
  lines(x=X, y=emaTA(Y), col="black") }
plot(tS3, plot.type="m", col = .colorwheelPalette(3), panel=lines2, 
  grid=TRUE, at="pretty")


###################################################
### code chunk number 67: margins
###################################################
# Plot:
# - oma stands for 'Outer Margin Area'
# - mar represents the 'figure Margins'  
# - The default size is c(5,4,4,2) + 0.1
# - The axes tick marks will go in the first lines 
par(mfrow=c(1, 1))
par(oma=c(3,3,3,3))  # all sides have 3 lines of space  
par(mar=c(5,4,4,2) + 0.1)    
plot(x=1:10, y=1:10, type="n", xlab="X", ylab="Y")   
  
# Add Text tot the Plot Part - red
text(5,5, "Plot", col="red", cex=2)  
text(5,4, "text(5,5, \"Plot\", col=\"red\", cex=2)", col="red", cex=1)  
box("plot", col="red", lwd=2) 

# Add text to thebThe Figure Part - grey
mtext("Margins", side=3, line=2, cex=1.5, col="grey")  
mtext("par(mar=c(5,4,4,2) + 0.1)", side=3, line=1, cex=1, col="grey")  
mtext("Line 0", side=3, line=0, adj=1.0, cex=1, col="grey")  
mtext("     1", side=3, line=1, adj=1.0, cex=1, col="grey")  
mtext("Line 2", side=3, line=2, adj=1.0, cex=1, col="grey")  
mtext("Line 3", side=3, line=3, adj=1.0, cex=1, col="grey")  
mtext("Line 0", side=2, line=0, adj=1.0, cex=1, col="grey")  
mtext("Line 1", side=2, line=1, adj=1.0, cex=1, col="grey")  
mtext("Line 2", side=2, line=2, adj=1.0, cex=1, col="grey")  
mtext("Line 3", side=2, line=3, adj=1.0, cex=1, col="grey")  
box("figure", col="grey")  

# The title will fit in the third line on the top of the graph.   
title("Ttitle - Third Line") 
   
# Note 'outer=TRUE' moves us from the figure to the outer margins.  
mtext("Outer Margin Area", side=1, line=1, cex=1.8, col="brown", outer=TRUE)  
mtext("par(oma=c(3,3,3,3))", side=1, line=2, cex=1, col="orange", outer=TRUE)  
mtext("Line 0", side=1, line=0, adj=0.0, cex=0.8, col="orange", outer=TRUE)  
mtext("Line 1", side=1, line=1, adj=0.0, cex=1, col="orange", outer=TRUE)  
mtext("Line 2", side=1, line=2, adj=0.0, cex=1, col="orange", outer=TRUE)  
box("outer", col="orange") 


###################################################
### code chunk number 68: marginsPlot
###################################################
# Plot:
# - oma stands for 'Outer Margin Area'
# - mar represents the 'figure Margins'  
# - The default size is c(5,4,4,2) + 0.1
# - The axes tick marks will go in the first lines 
par(mfrow=c(1, 1))
par(oma=c(3,3,3,3))  # all sides have 3 lines of space  
par(mar=c(5,4,4,2) + 0.1)    
plot(x=1:10, y=1:10, type="n", xlab="X", ylab="Y")   
  
# Add Text tot the Plot Part - red
text(5,5, "Plot", col="red", cex=2)  
text(5,4, "text(5,5, \"Plot\", col=\"red\", cex=2)", col="red", cex=1)  
box("plot", col="red", lwd=2) 

# Add text to thebThe Figure Part - grey
mtext("Margins", side=3, line=2, cex=1.5, col="grey")  
mtext("par(mar=c(5,4,4,2) + 0.1)", side=3, line=1, cex=1, col="grey")  
mtext("Line 0", side=3, line=0, adj=1.0, cex=1, col="grey")  
mtext("     1", side=3, line=1, adj=1.0, cex=1, col="grey")  
mtext("Line 2", side=3, line=2, adj=1.0, cex=1, col="grey")  
mtext("Line 3", side=3, line=3, adj=1.0, cex=1, col="grey")  
mtext("Line 0", side=2, line=0, adj=1.0, cex=1, col="grey")  
mtext("Line 1", side=2, line=1, adj=1.0, cex=1, col="grey")  
mtext("Line 2", side=2, line=2, adj=1.0, cex=1, col="grey")  
mtext("Line 3", side=2, line=3, adj=1.0, cex=1, col="grey")  
box("figure", col="grey")  

# The title will fit in the third line on the top of the graph.   
title("Ttitle - Third Line") 
   
# Note 'outer=TRUE' moves us from the figure to the outer margins.  
mtext("Outer Margin Area", side=1, line=1, cex=1.8, col="brown", outer=TRUE)  
mtext("par(oma=c(3,3,3,3))", side=1, line=2, cex=1, col="orange", outer=TRUE)  
mtext("Line 0", side=1, line=0, adj=0.0, cex=0.8, col="orange", outer=TRUE)  
mtext("Line 1", side=1, line=1, adj=0.0, cex=1, col="orange", outer=TRUE)  
mtext("Line 2", side=1, line=2, adj=0.0, cex=1, col="orange", outer=TRUE)  
box("outer", col="orange") 


###################################################
### code chunk number 69: prettyAppendix
###################################################
FORMAT <- tS1@format
FORMAT
POSITIONS <- pretty(tS1)
POSITIONS
LABELS <- pretty(tS1)
LABELS


###################################################
### code chunk number 70: axTicks
###################################################
axTicksByTime <-
function (x, ticks.on = "auto", k = 1, labels = TRUE, format.labels = TRUE, 
    ends = TRUE, gt = 2, lt = 30) 
{
    if (timeBased(x)) x <- xts(rep(1, length(x)), x)
    tick.opts <- c("years", "months", "weeks", "days", "hours", "minutes", "seconds")
    tick.k.opts <- c(10, 5, 2, 1, 6, 1, 1, 1, 4, 2, 1, 30, 15, 1, 1)
    if (ticks.on %in% tick.opts) {
        cl <- ticks.on[1]
        ck <- k
    } else {
        tick.opts <- paste(rep(tick.opts, c(4, 2, 1, 1, 3, 3, 1)), tick.k.opts)
        is <- structure(rep(0, length(tick.opts)), .Names = tick.opts)
        for (i in 1:length(tick.opts)) 
        {
            y <- strsplit(tick.opts[i], " ")[[1]]
            ep <- endpoints(x, y[1], as.numeric(y[2]))
            is[i] <- length(ep) - 1
            if (is[i] > lt) break
        }
        nms <- rev(names(is)[which(is > gt & is < lt)])[1]
        cl <- strsplit(nms, " ")[[1]][1]
        ck <- as.numeric(strsplit(nms, " ")[[1]][2])
    }
    if (is.null(cl)) ep <- NULL else ep <- endpoints(x, cl, ck)
    if (ends) ep <- ep + c(rep(1, length(ep) - 1), 0)
    if (labels) 
    {
        if (is.logical(format.labels) || is.character(format.labels)) 
        {
            unix <- ifelse(.Platform$OS.type == "unix", TRUE, FALSE)
            time.scale <- periodicity(x)$scale
            fmt <- ifelse(unix, "%n%b%n%Y", "%b %Y")
            if (time.scale == "weekly" | time.scale == "daily") 
                fmt <- ifelse(unix, "%b %d%n%Y", "%b %d %Y")
            if (time.scale == "minute" | time.scale == "hourly") 
                fmt <- ifelse(unix, "%b %d%n%H:%M", "%b %d %H:%M")
            if (time.scale == "seconds") 
                fmt <- ifelse(unix, "%b %d%n%H:%M:%S", "%b %d %H:%M:%S")
            if (is.character(format.labels)) 
                fmt <- format.labels
            names(ep) <- format(index(x)[ep], fmt)
        } else {
           names(ep) <- as.character(index(x)[ep])
        }
      ep
   }
}


###################################################
### code chunk number 71: axTicks2
###################################################
ticks <- axTicksByTime(as.xts(tS1))
ticks


