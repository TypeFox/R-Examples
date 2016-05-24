


 x = tS1
 FinCenter = NULL
 type = NULL
 plot.type = c("multiple", "single")
 format = "auto"
 at = c("pretty", "chic")
 main <- xlab <- ylab <- ""; nm = colnames(x); log = ""
 col = 1; pch = 19; cex = 1; lty = 1; lwd = 1 
 grid = TRUE; frame.plot = TRUE
 xlim = NULL; ylim = NULL
 axes = TRUE; ann = TRUE; cex.axis = 1; cex.lab =1;
 yax.flip = FALSE
 mar.multi = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1)
 oma.multi = c(7.75, 1.1, 6.1, 1.1)
 
 






# Plot Function Extensions written by Diethel Wuertz
#   ... first Version 2014-05-12


###############################################################################
# 1    Standard Plots
# 1.1    Single Plots
# 1.2    Multiple Plots
# 1.2    Scatter Plots
# 2    Time Axis Layout
# 2.1    Pretty Axis Layout 
# 2.2    Chic Axis Layout
# 2.3    Tailored Axis Layout
# 3    Annotations
# 3.1    Adding Title and Labels
# 3.2    Removing Annotations
# 3.3    Changing Font Size
# 3.4    Flipping Value Axes
# 4    Decorations
# 4.1    Modifying Types
# 4.2    Changing Colors
# 4.3    Changing Line Styles
# 4.4    Changing Plot Symbols
# 4.5    Modifying Line Widths
# 4.6    Modifying Plot Symbol Sizes
###############################################################################

# First let us see what plot.ts can do in the multiple plot mode:

require(timeSeries)
tS1 <- 100 * cumulated(LPP2005REC[, 2])
tS2 <- 100 * cumulated(LPP2005REC[, 2:3])
tS3 <- 100 * cumulated(LPP2005REC[, 1:3])
tS6 <- 100 * cumulated(LPP2005REC[, 1:6])
tS7 <- 100 * cumulated(LPP2005REC[, 1:7])


# -----------------------------------------------------------------------------


                            

 
# 1.3 Scatter Plots:
 
mat <- getDataPart(tS2)
 
par(mfrow=c(2,2))
plot(mat[, 1], mat[, 2])  
plot(mat[, 1], mat[, 2], pch=19, cex=0.2)  
 

################################################################################
# 2. Time Axis Layout:
 
 
# Changing Time-Axis Size:


# One Column Multiple Plots - Each curve in its own Graph:
par(mfrow=c(1, 1))                
plot(tS3, at="chic", plot.type="m", cex.axis=0.8)    
 
par(mfrow=c(1, 1)) 
plot(tS3, at="chic", plot.type="m", cex.axis=1.1)   
 
 
# Two Columns Multiple Plots - Each curve in its own Graph:
 
par(mfrow=c(1, 1)) 
plot(tS6, at="chic", plot.type="m", cex.axis=0.8)  

par(mfrow=c(1, 1)) 
plot(tS6, at="chic", plot.type="m", cex.axis=1.1)  


################################################################################
# 3  Annotations
 
 
# ------------------------------------------------------------------------------
# 3.1  Adding Title and Labels


# Single Plot - All Curves in one Graph:   
par(mfrow=c(2, 2))
plot(tS1); title(main = "Index") 
plot(tS3, plot.type="s"); title(main = "Index") 
plot(tS3, plot.type="s"); title(main = "Index", xlab = "Date")                              
plot(tS6, plot.type="s"); title(main = "Index", xlab = "Date")  


# One Column Multiple Plots - Each curve in its own Graph:
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m"); title(main = "Index", xlab = "Date") 
 

# Two Column Multiple Plots:
par(mfrow=c(1, 1))
plot(tS6, plot.type="m"); title(main = "Index", xlab = "Date")  
 
 
# One Column Multiple Plots - User designed Title
par(mfrow=c(1, 1))                
plot(tS3, plot.type="m", at = "chic")    
mtext("Swiss Market", side=3, line=1, adj=-0.025)

 
# Two Column Multiple Plots - User designed Title
par(mfrow=c(1, 1)) 
plot(tS6, plot.type="m", at = "chic")  
mtext("Swiss Market", side=3, line=1, adj=-0.3)
mtext("Foreign Market", side=3, line=1, adj=0.7)

 
# ------------------------------------------------------------------------------
# 3.2 Remove all Annotations:
 
 
# Two Column Multiple Plots - Each curve in its own Graph:
par(mfrow=c(1, 1))
plot(tS6, plot.type="m", ann=FALSE)  
 

# Single Plot - All Curves in one Graph:
par(mfrow=c(2, 1), mar = c(4, 4, 1, 2) + 0.1)   
plot(tS1, at="chic", ann=FALSE) 
title(ylab = colnames(tS1), cex.axis=0.8, cex.lab=0.8)
plot(tS1, at="chic", ann=FALSE) 
title(ylab = colnames(tS1), cex.axis=1.2, cex.lab=1.2)


################################################################################
# 4  Decorations

 
# ------------------------------------------------------------------------------
# 4.1  Modifying Types

 
# "type"
par(mfrow=c(1,1))
plot(tS3, type = c("l", "p", "h"), plot.type="m")          
 
 
# ------------------------------------------------------------------------------
# 4.2  Changing Colors
 
 
# Selecting Colors:
par(mfrow=c(2, 2))
plot(tS3, col = 1, plot.type="s") 
plot(tS3, col = 1:3, plot.type="s") 
plot(tS3, col = c("blue", "orange", "brown"), plot.type="s")  
 
 
# ------------------------------------------------------------------------------ 
# 4.3  Changing Line Styles
 
 
# Single Plot:
par(mfrow=c(1,1))
plot(tS3, lty = 3:1, plot.type="s")  
 
 
# One Column Multiple Plot:
par(mfrow=c(1,1))
plot(tS3, lty = 3:1, plot.type="m")
 
 
# ------------------------------------------------------------------------------
# 4.4 Changing Plot Symbols

 
par(mfrow=c(1,1))
plot(tS3, pch = c(17, 18, 19), plot.type="m")  
 
 
# 4.5 Modifying Line Widths
 
 
par(mfrow=c(1,1))
plot(tS3, lwd = c(17, 18, 19), plot.type="m")  



par(mfrow=c(2, 2))
plot(tS3, type=rep("p", 3), cex = rep(1.2, 3), 
  plot.type="s")
plot(tS3, type = rep("p", 3), pch = 1:3, col = c("blue", "orange", "brown"), 
  plot.type="s") 
plot(tS3, type = "p", pch = 19, col = c("blue", "orange", "brown"), 
  plot.type="s")                            



# -----------------------------------------------------------------------------

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 chart.TimeSeries
function (R, auto.grid = TRUE, xaxis = TRUE, yaxis = TRUE, yaxis.right = FALSE, 
    type = "l", lty = 1, lwd = 2, main = NULL, ylab = NULL, xlab = "Date", 
    date.format.in = "%Y-%m-%d", date.format = NULL, xlim = NULL, 
    ylim = NULL, element.color = "darkgray", event.lines = NULL, 
    event.labels = NULL, period.areas = NULL, event.color = "darkgray", 
    period.color = "aliceblue", colorset = (1:12), pch = (1:12), 
    legend.loc = NULL, ylog = FALSE, cex.axis = 0.8, cex.legend = 0.8, 
    cex.lab = 1, cex.labels = 0.8, cex.main = 1, major.ticks = "auto", 
    minor.ticks = TRUE, grid.color = "lightgray", grid.lty = "dotted", 
    xaxis.labels = NULL, ...) 
{
    y = checkData(R)
    columns = ncol(y)
    rows = nrow(y)
    columnnames = colnames(y)
    if (is.null(date.format)) {
        freq = periodicity(y)
        yr_eq <- ifelse(format(index(first(y)), format = "%Y") == 
            format(index(last(y)), format = "%Y"), TRUE, FALSE)
        switch(freq$scale, seconds = {
            date.format = "%H:%M"
        }, minute = {
            date.format = "%H:%M"
        }, hourly = {
            date.format = "%d %H"
        }, daily = {
            if (yr_eq) date.format = "%b %d" else date.format = "%Y-%m-%d"
        }, weekly = {
            if (yr_eq) date.format = "%b %d" else date.format = "%Y-%m-%d"
        }, monthly = {
            if (yr_eq) date.format = "%b" else date.format = "%b %y"
        }, quarterly = {
            if (yr_eq) date.format = "%b" else date.format = "%b %y"
        }, yearly = {
            date.format = "%Y"
        })
    }
    rownames = as.Date(xts:::time.xts(y))
    rownames = format(strptime(rownames, format = date.format.in), 
        date.format)
    time.scale = periodicity(y)$scale
    ep = axTicksByTime(y, major.ticks, format.labels = date.format)
    logaxis = ""
    if (ylog) {
        logaxis = "y"
    }
    plot.new()
    if (is.null(xlim[1])) 
        xlim = c(1, rows)
    if (is.null(ylim[1])) {
        ylim = as.numeric(range(y, na.rm = TRUE))
    }
    plot.window(xlim, ylim, xaxs = "r", log = logaxis)
    if (is.null(ylab)) {
        if (ylog) 
            ylab = "ln(Value)"
        else ylab = "Value"
    }
    if (ylog) 
        dimensions = 10^par("usr")
    else dimensions = par("usr")
    if (!is.null(period.areas)) {
        period.dat = lapply(period.areas, function(x, y) c(first(index(y[x])), 
            last(index(y[x]))), y = y)
        period.ind = NULL
        for (period in 1:length(period.dat)) {
            if (!is.na(period.dat[[period]][1])) {
                period.ind = list(grep(period.dat[[period]][1], 
                  index(y)), grep(period.dat[[period]][2], index(y)))
                rect(period.ind[1], dimensions[3], period.ind[2], 
                  dimensions[4], col = period.color, border = NA)
            }
        }
    }
    if (auto.grid) {
        abline(v = ep, col = grid.color, lty = grid.lty)
        grid(NA, NULL, col = grid.color)
    }
    abline(h = 0, col = element.color)
    if (!is.null(event.lines)) {
        event.ind = NULL
        for (event in 1:length(event.lines)) {
            event.ind = c(event.ind, grep(event.lines[event], 
                rownames))
        }
        number.event.labels = ((length(event.labels) - length(event.ind) + 
            1):length(event.labels))
        abline(v = event.ind, col = event.color, lty = 2)
        if (!is.null(event.labels)) {
            text(x = event.ind, y = ylim[2], label = event.labels[number.event.labels], 
                offset = 0.2, pos = 2, cex = cex.labels, srt = 90, 
                col = event.color)
        }
    }
    if (length(lwd) < columns) 
        lwd = rep(lwd, columns)
    if (length(lty) < columns) 
        lty = rep(lty, columns)
    if (length(pch) < columns) 
        pch = rep(pch, columns)
    for (column in columns:1) {
        lines(1:rows, y[, column], col = colorset[column], lwd = lwd[column], 
            pch = pch[column], lty = lty[column], type = type, 
            ...)
    }
    if (xaxis) {
        if (minor.ticks) 
            axis(1, at = 1:NROW(y), labels = FALSE, col = "#BBBBBB")
        label.height = cex.axis * (0.5 + apply(t(names(ep)), 
            1, function(X) max(strheight(X, units = "in")/par("cin")[2])))
        if (is.null(xaxis.labels)) 
            xaxis.labels = names(ep)
        else ep = 1:length(xaxis.labels)
        axis(1, at = ep, labels = xaxis.labels, las = 1, lwd = 1, 
            mgp = c(3, label.height, 0), cex.axis = cex.axis)
        title(xlab = xlab, cex = cex.lab)
    }
    if (yaxis) 
        if (yaxis.right) 
            axis(4, cex.axis = cex.axis, col = element.color, 
                ylog = ylog)
        else axis(2, cex.axis = cex.axis, col = element.color, 
            ylog = ylog)
    box(col = element.color)
    if (!is.null(legend.loc)) {
        legend(legend.loc, inset = 0.02, text.col = colorset, 
            col = colorset, cex = cex.legend, border.col = element.color, 
            lty = lty, lwd = 2, bg = "white", legend = columnnames, 
            pch = pch)
    }
    if (is.null(main)) 
        main = columnnames[1]
    title(ylab = ylab, cex = cex.lab)
    title(main = main, cex = cex.main)
}
 
 
 

