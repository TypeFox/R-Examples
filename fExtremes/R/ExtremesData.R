
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
# FUNCTION             EXPLORATIVE DATA ANALYSIS:
#  emdPlot              Creates an empirical distribution plot
#  qqparetoPlot         Creates exploratory QQ plot for EV analysis
#  mePlot               Creates a sample mean excess plot
#   mxfPlot             Creates another view of a sample mean excess plot
#   mrlPlot             Returns a mean residual life plot with confidence levels
#  recordsPlot          Plots records development
#   ssrecordsPlot       Plots records development of data subsamples
#  msratioPlot          Plots ratio of maximums and sums
#  sllnPlot             Verifies Kolmogorov's Strong Law of large numbers
#  lilPlot              Verifies Hartman-Wintner's Law of the iterated logarithm
#  xacfPlot             Plots autocorrelations of exceedences
################################################################################


emdPlot = 
function(x, doplot = TRUE, plottype = c("xy", "x", "y", " "),
labels = TRUE, ...)
{   # A function imported from R-package evir

    # Description:
    #   Plots empirical distribution function
    
    # Arguments:
    #   x - any object which can be transformed by the function
    #       as.vector() into a numeric vector 
    #   doplot - a logical flag, should a pot be returned ?
    #   plottype - which axes should be on a log scale: 
    #       "x" denotes x-axis only; "y" denotes y-axis only,
    #       "xy" || "yx" both axes, "" denotes neither of the 
    #       axis

    # FUNCTION:
        
    # Convert Type:
    x = as.vector(x)
    
    # Settings:
    plottype = match.arg(plottype)
         
    # Convert x to a vector, if the input is a data.frame.
    if (is.data.frame(x)) x = x[, 1] 
    xs = x = sort(as.numeric(x))
    ys = y = 1 - ppoints(x)
    
    if (plottype == "x") {
        xs = x[x > 0]
        ys = y[x > 0] 
    }
    if (plottype == "y") {
        xs = x[y > 0]
        ys = y[y > 0] 
    }
    if (plottype == "xy") {
        xs = x[x > 0 & y > 0]
        ys = y[x > 0 & y > 0] 
    }
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "x"
            ylab = "1-F(x)"
            main = "Empirical Distribution" 
            if (plottype == "xy") main = paste("log-log", main)
            if (plottype ==  "x") main = paste("log-lin", main)
            if (plottype ==  "y") main = paste("lin-log", main)
            if (plottype ==   "") main = paste("lin-lin", main)
        } else {
            xlab = ""
            ylab = ""
            main = "" 
        }   
        if (labels) {
            plot(xs, ys, pch = 19, col = "steelblue",
                log = plottype, xlab = xlab, ylab = ylab, main = main, ...) 
            grid()
        } else {
            plot(xs, ys, 
                log = plottype, xlab = xlab, ylab = ylab, main = main, ...) 
        }
    }           
    
    # Result:
    result = data.frame(x, y)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)   
}


# ------------------------------------------------------------------------------


qqparetoPlot = 
function(x, xi = 0, trim = NULL, threshold = NULL, doplot = TRUE, 
labels = TRUE, ...)
{   # A function imported from R-package evir
    
    # Description:
    #   Creates an exploratory QQ-plot for Extreme Value Analysis.

    # Arguments:
    #   x - any object which can be transformed by the function
    #       as.vector() into a numeric vector 
    #   doplot - a logical flag, should a plot be returned ?
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Convert x to a vector, if the input is a data.frame.
    if(is.data.frame(x)) x = x[, 1] 
    
    # qPlot:
    x = as.numeric(x)
    if (!is.null(threshold)) x = x[x >= threshold]
    if (!is.null(trim)) x = x[x < trim]
    if (xi == 0) {
        y = qexp(ppoints(x)) 
    }
    if( xi != 0) {
        y = qgpd(ppoints(x), xi = xi) 
    }
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "Ordered Data"
            ylab = "Quantiles"
            if (xi == 0) {
                ylab = paste("Exponential", ylab) 
            }
            if (xi != 0) {
                ylab = paste("GPD(xi=", xi, ") ",  ylab, sep = "") 
            }
            main = "Exploratory QQ Plot" 
        } else {
            xlab = ""
            ylab = ""
            main = "" 
        }
        z = sort(x)
        plot(z, y, pch = 19, col = "steelblue", xlab = xlab, 
            ylab = ylab, main = main, ...)
        rug(z, ticksize = 0.01, side = 3)
        rug(y, ticksize = 0.01, side = 4)
        abline(lsfit(z, y)) 
        if (labels) {
            grid()
            text = paste("xi =", as.character(round(xi, 3))) 
            mtext(text, side = 4, adj = 0, cex = 0.7)
        }
    }
    
    # Result:
    result = data.frame(x = sort(x), y)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


# ------------------------------------------------------------------------------


mxfPlot = 
function (x, u = quantile(x, 0.05), doplot = TRUE, labels = TRUE, ...)     
{   # A function written by Diethelm Wuertz
    
    # Description:
    #   Creates a simple mean excess function plot.
    
    # Arguments:
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Convert x to a vector, if the input is a data.frame.
    if(is.data.frame(x)) x = x[, 1] 
    
    # mxf:
    tail = length(x[x < u])/length(x)
    u = rev(sort(x))
    n = length(x)
    u = u[1:floor(tail*n)]
    n = length(u)
    e = (cumsum(u)-(1:n)*u)/(1:n)
    
    # Plot
    if (doplot) {
        if (labels) {
            xlab = "Threshold: u"
            ylab = "Mean Excess: e"
            main = "Mean Excess Function" 
        } else {
            main = xlab = ylab = ""
        }
        plot (u, e, pch = 19, col = "steelblue",  
             xlab = xlab, ylab = ylab, main = main, ...) 
        if (labels) grid()
    }
    
    # Result:
    result = data.frame(threshold = u, excess = e)
    
    # Return Values:
    if (doplot) return(invisible(result)) else return(result)
}


# ------------------------------------------------------------------------------


mrlPlot = 
function(x, ci = 0.95, umin = mean(x), umax = max(x), nint = 100, 
doplot = TRUE, plottype = c("autoscale", ""), labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Create a mean residual life plot with
    #   confidence intervals.
    
    # Arguments:
    
    # References:
    #   A function originally written by S. Coles
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Settings:
    plottype = plottype[1]
    
    # Convert x to a vector, if the input is a data.frame.
    if (is.data.frame(x)) x = x[,1] 
    sx = xu = xl = rep(NA, nint)
    u = seq(umin, umax, length = nint)
    for (i in 1:nint) {
        x = x[x >= u[i]]
        sx[i] = mean(x - u[i])
        sdev = sqrt(var(x))
        n = length(x)
        xu[i] = sx[i] + (qnorm((1 + ci)/2) * sdev) / sqrt(n)
        xl[i] = sx[i] - (qnorm((1 + ci)/2) * sdev) / sqrt(n) 
    }
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "Threshold: u"
            ylab = "Mean Excess: e"
            main = "Mean Residual Live Plot" 
        } else {
            main = xlab = ylab = ""
        }
        if (plottype == "autoscale") {
            ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)]))
            plot(u, sx, type = "o", pch = 19, col = "steelblue",
                xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...) 
        } else {
            plot(u[!is.na(xl)], sx[!is.na(xl)], type = "o", 
                pch = 19, col = "steelblue",
                xlab = xlab, ylab = ylab, main = main, ...) 
        } 
        lines(u[!is.na(xl)], xl[!is.na(xl)], col = "brown")
        lines(u[!is.na(xu)], xu[!is.na(xu)], col = "brown")
        if (labels) {
            grid()
            text = paste("ci =", as.character(round(ci, 3))) 
            mtext(text, side = 4, adj = 0, cex = 0.7)
        } 
    }
    
    # Result
    result = data.frame(threshold = u, mrl = sx)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


# ------------------------------------------------------------------------------


mePlot = 
function(x, doplot = TRUE, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Create a Mean Excess Plot
    
    # Arguments:
    #   x - an univariate time series object or any other object which 
    #       can be transformed by the function as.vector() into a numeric
    #       vector.
    #   doplot - a logical flag, should a plot be created?
    #   labels - a logical flag, should the plot be automatically labeld?
    #       If TRUE, then default values to xlab, ylab, main, pch and col 
    #       are assigned.
    
    # Reference:
    #   A function imported from R-package evir
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Settings:
    omit = 0
    
    # Internal Function:
    myrank = function(x, na.last = TRUE){
        ranks = sort.list(sort.list(x, na.last = na.last))
        if (is.na(na.last))
            x = x[!is.na(x)]
        for (i in unique(x[duplicated(x)])) {
            which = x == i & !is.na(x)
            ranks[which] = max(ranks[which]) 
        }
        ranks 
    }
    
    # Convert x to a vector, if the input is a data.frame.
    if(is.data.frame(x)) x = x[, 1] 
    x = as.numeric(x)
    n = length(x)
    x = sort(x)
    n.excess = unique(floor(length(x) - myrank(x)))
    points = unique(x)
    nl = length(points)
    n.excess = n.excess[-nl]
    points = points[-nl]
    excess = cumsum(rev(x))[n.excess] - n.excess * points
    y = excess/n.excess
    xx = points[1:(nl-omit)] 
    yy = y[1:(nl-omit)]
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "Threshold: u"
            ylab = "Mean Excess: e"
            main = "Mean Excess Plot" 
            plot(xx, yy, pch = 19, col = "steelblue", 
                xlab = xlab, ylab = ylab, main = main, ...) 
            grid()
        } else {
            plot(xx, yy, ...) 
        }
    }
    
    # Results:
    result = data.frame(threshold = xx, me = yy)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)   
}


# -----------------------------------------------------------------------------


recordsPlot = 
function(x, ci = 0.95, doplot = TRUE, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a records plot.
    
    # Note:
    #   A function imported from R-package evir,
    #   original name in EVIR: records

    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Settings:
    conf.level = ci
    
    # Convert x to a vector, if the input is a data.frame.
    if (is.data.frame(x)) x = x[,1] 
    
    # Records:
    record = cummax(x)
    expected = cumsum(1/(1:length(x)))
    se = sqrt(expected - cumsum(1/((1:length(x))^2)))
    trial = (1:length(x))[!duplicated(record)]
    record = unique(record)
    number = 1:length(record)
    expected = expected[trial]
    se = se[trial]
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "Trials"
            ylab = "Records"
            main = "Plot of Record Development" 
        } else {
            xlab = ""
            ylab = ""
            main = "" 
        }     
        ci = qnorm(0.5 + conf.level/2)
        upper = expected + ci * se
        lower = expected - ci * se
        lower[lower < 1] = 1
        yr = range(upper, lower, number)    
        plot(trial, number, log = "x", ylim = yr, pch = 19,
            col = "steelblue", xlab = xlab, ylab = ylab, 
            main = main, ...) 
        lines(trial, expected)
        lines(trial, upper, lty = 2, col = "brown")
        lines(trial, lower, lty = 2, col = "brown") 
        if (labels) {
            grid()
            text = paste("ci =", as.character(conf.level)) 
            mtext(text, side = 4, adj = 0, col = "grey", cex = 0.7)
        } 
    }
        
    # Result:
    result = data.frame(number, record, trial, expected, se)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


# ------------------------------------------------------------------------------


ssrecordsPlot = 
function (x, subsamples = 10, doplot = TRUE, plottype = c("lin", "log"),
labels = TRUE,  ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a plot of records on subsamples.
    
    # note:
    #   Changes:
    #   2003/09/06 - argument list made consistent
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Plot type:
    plottype = match.arg(plottype)
    
    # Labels:
    xlab = ylab = main = ""    
    
    # Records:
    save = x 
    cluster = floor(length(save)/subsamples)
    records = c()
    for (i in 1:subsamples) {
        x = save[((i-1)*cluster+1):(i*cluster)]
        y = 1:length(x)
        u = x[1]
        v = x.records = 1
        while (!is.na(v)) {
            u = x[x > u][1]
            v = y[x > u][1]
            if(!is.na(v)) x.records = c(x.records, v) 
        }   
        if (i == 1) {
            nc = 1:length(x)
            csmean = cumsum(1/nc)
            cssd = sqrt(cumsum(1/nc-1/(nc*nc)))
            ymax = csmean[length(x)] + 2*cssd[length(x)]
            # Plot:
            if (doplot) {
                if (plottype == "log") {
                    nc = log(nc)
                }
                if (labels) {
                    if (plottype == "lin") xlab = "n"
                    if (plottype == "log") xlab = "log(n)"
                    ylab = "N(n)" 
                    main = "Subsample Records Plot"
                }
                plot (nc, csmean+cssd, type = "l", ylim = c(0, ymax),
                    lty = 2, col = "brown", xlab = xlab, ylab = ylab, 
                    main = main, ...) 
                lines(nc, csmean)  
                lines(nc, csmean-cssd, lty = 2, col = "brown") 
                if (labels) {
                    grid() 
                    text = paste("subsamples =", as.character(subsamples)) 
                    mtext(text, side = 4, adj = 0, col = "grey", cex = 0.7)
                }
            }
        } 
        y.records = 1:length(x.records)
        x.records = x.records[y.records < ymax]
        if (doplot) {
            if (plottype == "log") {
                x.records = log(x.records)
            }
            points(x.records, y.records[y.records<ymax], 
               pch = i, col = "steelblue") 
        }
        records[i] = y.records[length(y.records)]
    }
    
    # Result:
    subsample = 1:subsamples
    result = data.frame(subsample, records)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


# ------------------------------------------------------------------------------


msratioPlot = 
function (x, p = 1:4, doplot = TRUE, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a Plot of maximum and sum ratio.
    
    # Arguments:
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # Settings:
    P = length(p)
    plottype = "autoscale"
    
    # Convert x to a vector, if the input is a data.frame.
    if(is.data.frame(x)) x = x[, 1] 
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "Trials"
            ylab = "Records"
            main = "Plot of Maximum/Sum Ratio" 
        } else {
            xlab = ""
            ylab = ""
            main = "" 
        }             
        if (plottype[1] == "autoscale") {
            ylim = c(0, 1)
            plot(c(0, length(x)), y = ylim, xlab = xlab, 
                    ylab = ylab, main = main, type = "n", ...) 
        } else {
            plot(c(0, length(x)), xlab = xlab, 
                    ylab = ylab, main = main, type = "n", ...) 
        }
        if (labels) grid()
    }
    
    # Suppress warnings for points outside the frame:
    ratios = matrix(rep(0, times = length(x)*P), byrow = TRUE, ncol = P)
    if (doplot) par(err = -1)

    # Loop over all exponents p:
    i = 1
    for (q in p) {
        rnp = cummax(abs(x)^q) / cumsum(abs(x)^q)
        ratios[, i] = rnp
        i = i + 1
        if (doplot) lines (rnp, col = i, lty = i) 
    }

    # Result:
    result = data.frame(ratios)
    
    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


# ------------------------------------------------------------------------------


sllnPlot =  
function(x, doplot = TRUE, labels = TRUE, ...)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Verifies Kolmogorov's strong law of large numbers
    
    # Arguments:
    #   x - sequence of iid non-degenerate rvs.
    
    # References:
    #   Embrechts et al. p. 61, Theorem 2.1.3
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # SLLN:
    if (is.null(mean)) mean = mean(cumsum(x)/(1:length(x)))
    nx  =  length(x)
    
    # Plot:
    y  = cumsum(x)/(1:nx)
    mean = mean(x)
    if (doplot) {
        if (labels) {
            xlab = "n"
            ylab = "x"
            main = "SLLN" 
        } else {
            xlab = ""
            ylab = ""
            main = "" 
        }             
        plot(y, xlab = xlab, ylab = ylab, type = "l", main = main, 
            col = "steelblue", ...)
        lines(c(0, nx), c(mean, mean), col = "brown")
        if (labels) grid()
    }
    
    # Return Value:
    invisible(y)
}


# ------------------------------------------------------------------------------


lilPlot =  
function (x, doplot = TRUE, labels = TRUE, ...)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Verifies Hartman-Wintner's Law of the iterated logarithm
            
    # Arguments:
    #   x - sequence of iid non-degenerate rvs.

    # References:
    #   Embrechts et al. p. 67. Theorem 2.1.13
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    
    # LIL:
    lx = length(x)
    nx = 3:lx
    fact = sqrt(2*nx*log(log(nx)))
    mu = mean(x)  
    sdev = sqrt(var(x))
    y = (cumsum(x)[-(1:2)]-mu*nx)/(fact*sdev)
    
    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "n"
            ylab = "x"
            main = "LIL" 
        } else {
            xlab = ""
            ylab = ""
            main = "" 
        }        
        plot(x = nx, y = y, xlab = "n", ylab = "x", 
            ylim = range(y[!is.na(y)], -1, 1), type = "l", 
            main = main, col = "steelblue", ...)
        lines(c(0,lx), c(1,1), col = "brown")
        lines(c(0,lx), c(-1,-1), col = "brown")
        if (labels) grid()
    }
    
    # Return Value:
    invisible(y)
}


# ------------------------------------------------------------------------------


xacfPlot = 
function(x, u = quantile(x, 0.95), lag.max = 15, doplot = TRUE, 
which = c("all", 1, 2, 3, 4), labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates plots of exceedences, one for the
    #   heights and one for the distances.
    
    # Arguments:
    
    # FUNCTION:
    
    # Convert Type:
    x = as.vector(x)
    # which = match.arg(which)
    which = as.character(which[1])
    stopifnot(which == "all" | 
        which == "1" | which == "2" | which == "3" | which == "4")
    which = as.character(which)
    
    # Settings:
    if (labels) {
        xlab = c("Index", "Lag")
        ylab = c("Heights", "Distances", "ACF")
        main = c("Heights over Threshold", "Distances between Heights", 
            "Series Heights", "Series Distances") 
    } else {
        xlab = c("", "")
        ylab = c("", "", "")
        main = c("", "", "", "")
    } 
    
    # Heights/Distances:
    Heights = (x-u)[x > u]
    Distances = diff((1:length(x))[x > u])
    
    # Plot:
    if (doplot) {
        if (which == "all" | which == "1")
        plot (Heights, type = "h", xlab = xlab[1], ylab = ylab[1], 
            main = main[1], ...)
        if (which == "all" | which == "2")
        plot (Distances, type = "h", xlab = xlab[1], ylab = ylab[2], 
            main = main[2], ...) 
    }
    
    # Correlations:
    if (which == "all" | which == "3")
    Heights = as.vector(acf(Heights, lag.max=lag.max, plot = doplot, 
        xlab = xlab[2], ylab = ylab[3], main = main[3], ...)$acf)
    if (which == "all" | which == "4")
    Distances = as.vector(acf(Distances, lag.max=lag.max, plot = doplot, 
        xlab = xlab[2], ylab = ylab[3], main = main[4], ...)$acf)

    # Result:
    if (which == "all") {
        lag = as.vector(0:(lag.max))
        result = data.frame(lag, Heights, Distances) 
    } else {
        result = NULL
    }

    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


################################################################################

