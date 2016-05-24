
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
# FUNCTION:             MDA ESTIMATORS:
#  hillPlot              Plot Hill's estimator
#  shaparmPlot           Pickands, Hill & Decker-Einmahl-deHaan Estimator
#   shaparmPickands      Auxiliary function called by shaparmPlot
#   shaparmHill           ... called by shaparmPlot
#   shaparmDehaan         ... called by shaparmPlot
################################################################################


hillPlot =
function(x, start = 15, ci = 0.95,
doplot = TRUE, plottype = c("alpha", "xi"), labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots the results from the Hill Estimator.

    # Note:
    #   Code partly adapted from R package evir

    # Examples:
    #   par(mfrow = c(2, 2))
    #   hillPlot(gevSim(n=1000), plottype = "alpha")
    #   hillPlot(gevSim(n=1000), plottype = "xi")
    #   NYI: hillPlot(gevSim(n=1000), plottype = "alpha", reverse = TRUE)
    #   NYI: hillPlot(gevSim(n=1000), plottype = "xi", reverse = TRUE)
    #   hillPlot(gevSim(n=1000), plottype = "alpha", doplot = FALSE)
    #   hillPlot(gevSim(n=1000), plottype = "xi", doplot = FALSE)

    # Check Type:
    stopifnot(NCOL(x)==1)
    x = as.vector(x)

    # Settings:
    reverse = FALSE
    option = match.arg(plottype)
    data = x

    # MDA:
    ordered = rev(sort(data))
    ordered = ordered[ordered > 0]
    n = length(ordered)
    k = 1:n
    loggs = log(ordered)
    avesumlog = cumsum(loggs)/(1:n)
    xihat = c(NA, (avesumlog - loggs)[2:n])
    y = switch(option,
        alpha = 1/xihat,
        xi = xihat)
    ses = y / sqrt(k)
    x = trunc(seq(from = min(n, length(data)), to = start))
    y = y[x]
    qq <- qnorm(1 - (1 - ci)/2)
    u <- y + ses[x] * qq
    l <- y - ses[x] * qq
    yrange <- range(u, l)
    if (reverse) index = -x else index = x

    # Plot:
    if (doplot) {
        plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "",
            axes = FALSE, ...)
        pos = floor(seq(1, length(index), length = 10))
        axis(1, at = index[pos], labels = paste(x[pos]), tick = TRUE)
        axis(2)
        threshold = signif(findThreshold(data, x), 3)
        axis(3, at = index[pos], labels = paste(format(threshold[pos])),
            tick = TRUE)
        box()
        lines(index, u, lty = 2, col = "steelblue")
        lines(index, l, lty = 2, col = "steelblue")
        if (labels) {
            title(xlab = "Order Statistics", ylab = option)
            mtext("Threshold", side = 3, line = 3)
        }
    }

    # Result:
    ans = list(x = index, y = y)
    control = data.frame(plottype = option[1], start = start, ci = ci,
        reverse = FALSE, row.names = "control")
    attr(ans, "control") = control

    # Return Value:
    if (doplot) return(invisible(ans)) else ans
}


# ------------------------------------------------------------------------------


shaparmPlot =
function(x, p = 0.01*(1:10), xiRange = NULL, alphaRange = NULL,
doplot = TRUE, plottype = c("both", "upper"))
{   # A function written by Diethelm Wuertz

    # Description:
    #   Displays Pickands, Einmal-Decker-deHaan, and Hill estimators

    # Example:
    #   par(mfcol=c(3,2)); shaparmPlot(as.timeSeries(data(daxRet)))
    #   shaparmPlot(as.timeSeries(data(daxRet)), doplot = FALSE)
    #   shaparmPlot(as.timeSeries(data(daxRet)), 0.005*(1:20))

    # FUNCTION:

    # Settings:
    x = as.vector(x)
    tails = p
    if (is.null(xiRange)) xiRange = c(-0.5, 1.5)
    if (is.null(alphaRange)) alphaRange = c(0, 10)
    plottype = match.arg(plottype)
    if (plottype == "both") bothTails = TRUE else bothTails = FALSE

    # Median Plot:
    index = which.min(abs(tails-median(tails)))
    DOPLOT = rep(FALSE, length(tails))
    DOPLOT[index] = TRUE
    selected.tail = tails[index]
    if (!doplot) DOPLOT[index] = FALSE

    # Which estimator ?
    which = c(TRUE, TRUE, TRUE)

    # Settings:
    select.doplot = which
    ylim1 = xiRange
    ylim2 = alphaRange
    z = rep(mean(ylim2), length(tails))
    ylim1 = xiRange
    ylim2 = alphaRange

    # Estimates:
    p1 = p2 = h1 = h2 = d1 = d2 = m1 = m2 = rep(0, length(tails))
    for ( i in (1:length(tails)) ) {
        tail = tails[i]
        # Plotting Staff:
        if (select.doplot[1]) {
            xi = shaparmPickands(x, tail, ylim1, doplot = FALSE,
                plottype = plottype)
            p1[i] = xi$xi[1]
            p2[i] = xi$xi[3]
        }
        if (select.doplot[2]) {
            xi = shaparmHill(x, tail, ylim1, doplot = FALSE,
                plottype = plottype)
            h1[i] = xi$xi[1]
            h2[i] = xi$xi[3]
        }
        if (select.doplot[3]) {
            xi = shaparmDEHaan(x, tail, ylim1, doplot = FALSE,
                plottype = plottype)
            d1[i] = xi$xi[1]
            d2[i] = xi$xi[3]
        }
    }

    # Plot Pickands' Summary:
    if (select.doplot[1] & doplot) {
        plot (tails, z, type = "n", xlab = "tail depth", ylab = "alpha",
            ylim = ylim2, main = "Pickands Summary")
        grid()
        abline(v = selected.tail, lty = 3)
        y1 = 1/p1
        x1 = tails [y1 > ylim2[1] & y1 < ylim2[2]]
        y1 = y1[y1 > ylim2[1] & y1 < ylim2[2]]
        points (x1, y1, col = "steelblue")
        lines(x1, y1, col = "steelblue")
        if (bothTails) {
            y1 = 1/p2
            x1 = tails [y1 > ylim2[1] & y1 < ylim2[2]]
            y1 = y1 [y1 > ylim2[1] & y1 < ylim2[2]]
            points (x1, y1, col = "brown")
            lines(x1, y1, col = "brown")
        }
    }

    # Plot Hill Summary:
    if (select.doplot[2] & doplot) {
        plot (tails, z, type = "n", xlab = "tail depth", ylab = "alpha",
            ylim = ylim2, main = "Hill Summary")
        grid()
        abline(v = selected.tail, lty = 3)
        y1 = 1/h1
        x1 = tails [y1 > ylim2[1] & y1 < ylim2[2]]
        y1 = y1 [y1 > ylim2[1] & y1 < ylim2[2]]
        points (x1, y1, col = "steelblue")
        lines(x1, y1, col = "steelblue")
        if (bothTails) {
            y1 = 1/h2
            x1 = tails [y1 > ylim2[1] & y1 < ylim2[2]]
            y1 = y1 [y1 > ylim2[1] & y1 < ylim2[2]]
            points (x1, y1, col = "brown")
            lines(x1, y1, col = "brown")
        }
    }

    # Plot Deckers-Einmahl-deHaan Summary
    if (select.doplot[3] & doplot) {
        plot (tails, z, type = "n", xlab = "tail depth", ylab = "alpha",
            ylim = ylim2, main = "Deckers-Einmahl-deHaan Summary")
        grid()
        abline(v = selected.tail, lty = 3)
        y1 = 1/d1
        x1 = tails [y1>ylim2[1] & y1<ylim2[2]]
        y1 = y1 [y1>ylim2[1] & y1<ylim2[2]]
        points (x1, y1, col = "steelblue")
        lines(x1, y1, col = "steelblue")
        if (bothTails) {
            y1 = 1/d2
            x1 = tails [y1 > ylim2[1] & y1 < ylim2[2]]
            y1 = y1 [y1 > ylim2[1] & y1 < ylim2[2]]
            points (x1, y1, col = "brown")
            lines(x1, y1, col = "brown")
        }
    }

    # Plot Estimates:
    resultUpper = resultLower = NULL
    for ( i in (1:length(tails)) ) {
        tail = tails[i]
        # Plotting Staff:
        if (select.doplot[1]) {
            xi = shaparmPickands(x, tail, ylim1, doplot = DOPLOT[i],
                plottype = plottype)
            p1[i] = xi$xi[1]
            p2[i] = xi$xi[3]
        }
        if (select.doplot[2]) {
            xi = shaparmHill(x, tail, ylim1, doplot = DOPLOT[i],
                plottype = plottype)
            h1[i] = xi$xi[1]
            h2[i] = xi$xi[3]
        }
        if (select.doplot[3]) {
            xi = shaparmDEHaan(x, tail, ylim1, doplot = DOPLOT[i],
                plottype = plottype)
            d1[i] = xi$xi[1]
            d2[i] = xi$xi[3]
        }
        resultUpper = rbind(resultUpper, c(tails[i], p1[i], h1[i], d1[i]))
        if (bothTails)
            resultLower = rbind(resultLower, c(tails[i], p2[i], h2[i], d2[i]))
    }
    colnames(resultUpper) = c("Upper", "Pickands", "Hill", "DEHaan")
    resultUpper = data.frame(resultUpper)
    if (bothTails) {
        colnames(resultLower) = c("Lower", "Pickands", "Hill", "DEHaan")
        resultLower = data.frame(resultLower)
    }

    # Result:
    ans = list(Upper = resultUpper)
    if (bothTails) ans$Lower = resultLower

    # Return Value:
    if (doplot) return(invisible(ans)) else ans
}


# ------------------------------------------------------------------------------


shaparmPickands =
function(x, p = 0.05, xiRange = NULL,
doplot = TRUE, plottype = c("both", "upper"), labels = TRUE, ...)
{   # A function written by Diethelm Wuertz

    # FUNCTION:

    # Order Residuals:
    x = as.vector(x)
    tail = p
    if (is.null(xiRange)) xiRange = c(-0.5, 1.5)
    yrange = xiRange
    plottype = match.arg(plottype)
    if (plottype == "both") bothTails = TRUE else bothTails = FALSE
    ordered1 = rev(sort(abs(x[x < 0])))
    if (bothTails) ordered2 = rev(sort(abs(x[x > 0])))
    n1 = length(ordered1)
    if (bothTails) n2 = length(ordered2)
    ordered1 = ordered1[1:floor(tail*n1)]
    if (bothTails) ordered2 = ordered2[1:floor(tail*n2)]
    n1 = length(ordered1)
    if (bothTails) n2 = length(ordered2)

    # Pickands Estimate:
    k1 = 1:(n1%/%4)
    if (bothTails) k2 = 1:(n2%/%4)
    pickands1 = log ((c(ordered1[k1])-c(ordered1[2*k1])) /
        (c(ordered1[2*k1])-c(ordered1[4*k1]))) / log(2)
    if (bothTails) pickands2 = log ((c(ordered2[k2])-c(ordered2[2*k2])) /
        (c(ordered2[2*k2])-c(ordered2[4*k2]))) / log(2)

    # Prepare Plot:
    y1 = pickands1[pickands1 > yrange[1] & pickands1 < yrange[2]]
    x1 = log10(1:length(pickands1))[pickands1 > yrange[1] &
        pickands1 < yrange[2]]
    if (bothTails) {
        y2 = pickands2[pickands2 > yrange[1] & pickands2 < yrange[2]]
        x2 = log10(1:length(pickands2))[pickands2 > yrange[1] &
            pickands2 < yrange[2]]
    }
    # Labels:
    if (labels) {
        main = "Pickands Estimator"
        xlab = "log scale"
        ylab = "xi"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:
    if (doplot) {
        par(err = -1)
        plot (x1, y1, xlab = xlab, ylab = ylab, ylim = yrange,
            main = main, type = "n")
        title(sub = paste("tail depth:", as.character(tail)))
            lines(x1, y1, type = "p", pch = 2, col = "steelblue")
        if (bothTails) lines(x2, y2, type = "p", pch = 6, col = "brown")
        if (labels) grid()
    }

    # Calculate invers "xi":
    my1 = mean(y1, na.rm = TRUE)
    if (bothTails) my2 = mean(y2, na.rm = TRUE)
    sy1 = sqrt(var(y1, na.rm = TRUE))
    if (bothTails) sy2 = sqrt(var(y2, na.rm = TRUE))

    # Plot:
    if (doplot) {
        par(err = -1)
        lines(c(x1[1], x1[length(x1)]), c(my1,my1), type = "l",
            lty = 1, col = "steelblue")
        if (bothTails) lines(c(x2[1], x2[length(x2)]), c(my2, my2),
            type = "l", lty = 1, col = "brown")
    }

    # Result:
    result = list(xi = c(my1, sy1))
    if (bothTails) result = list(xi = c(my1, sy1, my2, sy2))

    # Return Result:
    result
}


# ------------------------------------------------------------------------------


shaparmHill =
function(x, p = 0.05, xiRange = NULL,
doplot = TRUE, plottype = c("both", "upper"), labels = TRUE, ...)
{   # A Function written by Diethelm Wuertz

    # ORDER RESIDUALS:
    x = as.vector(x)
    tail = p
    if (is.null(xiRange)) xiRange = c(-0.5, 1.5)
    yrange = xiRange
    plottype = match.arg(plottype)
    if (plottype == "both") bothTails = TRUE else bothTails = FALSE
    ordered1 = rev(sort(abs(x[x < 0])))
    if (bothTails) ordered2 = rev(sort(abs(x[x > 0])))
    n1 = length(ordered1)
    if (bothTails) n2 = length(ordered2)
    ordered1 = ordered1[1:floor(tail*n1)]
    if (bothTails) ordered2 = ordered2[1:floor(tail*n2)]
    n1 = length(ordered1)
    if (bothTails) n2 = length(ordered2)

    # HILLS ESTIMATE:
    hills1 = c((cumsum(log(ordered1))/(1:n1)-log(ordered1))[2:n1])
    if (bothTails) hills2 = c((cumsum(log(ordered2))/(1:n2) -
        log(ordered2))[2:n2])

    # PREPARE PLOT:
    y1 = hills1[hills1 > yrange[1] & hills1 < yrange[2]]
    x1 = log10(1:length(hills1))[hills1 > yrange[1] & hills1 < yrange[2]]
    if (bothTails) {
        y2 = hills2[hills2 > yrange[1] & hills2 < yrange[2]]
        x2 = log10(1:length(hills2))[hills2 > yrange[1] & hills2 < yrange[2]]
    }

    # Labels:
    if (labels) {
        main = "Hill Estimator"
        xlab = "log scale"
        ylab = "xi"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:
    if (doplot) {
        par(err = -1)
        plot (x1, y1, xlab = xlab, ylab = ylab, ylim = yrange,
            main = main, type="n")
        if (labels) title(sub = paste("tail depth:", as.character(tail)))
        lines(x1, y1, type = "p", pch = 2, col = "steelblue")
        if (bothTails) lines(x2, y2, type = "p", pch = 6, col = "brown")
        if (labels) grid()
    }

    # CALCULATE INVERSE XI:
    my1 = mean(y1, na.rm = TRUE)
    if (bothTails) my2 = mean(y2, na.rm = TRUE)
    sy1 = sqrt(var(y1, na.rm = TRUE))
    if (bothTails) sy2 = sqrt(var(y2, na.rm = TRUE))
    if (doplot) {
        par(err=-1)
        lines(c(x1[1], x1[length(x1)]), c(my1,my1), type="l",
            lty = 1, col = "steelblue")
        if (bothTails) lines(c(x2[1], x2[length(x2)]), c(my2,my2),
            type = "l",lty = 1, col = "brown")
    }

    # Result:
    result = list(xi = c(my1, sy1))
    if (bothTails) result = list(xi = c(my1, sy1, my2, sy2))

    # Return Result:
    result
}


# ------------------------------------------------------------------------------


shaparmDEHaan =
function(x, p = 0.05, xiRange = NULL,
doplot = TRUE, plottype = c("both", "upper"), labels = TRUE, ...)
{   # A function written by Diethelm Wuertz

    # ORDER RESIDUALS:
    x = as.vector(x)
    tail = p
    if (is.null(xiRange)) xiRange = c(-0.5, 1.5)
    yrange = xiRange
    plottype = match.arg(plottype)
    if (plottype == "both") bothTails = TRUE else bothTails = FALSE
    ordered1 = rev(sort(abs(x[x < 0])))
    if (bothTails) ordered2 = rev(sort(abs(x[x > 0])))
    n1 = length(ordered1)
    if (bothTails) n2 = length(ordered2)
    ordered1 = ordered1[1:floor(tail*n1)]
    if (bothTails) ordered2 = ordered2[1:floor(tail*n2)]
    n1 = length(ordered1)
    if (bothTails) n2 = length(ordered2)

    # DECKERS-EINMAHL-deHAAN ESTIMATE:
    ns0 = 1
    n1m = n1-1; ns1 = ns0; ns1p = ns1+1
    bod1 = c( cumsum(log(ordered1))[ns1:n1m]/(ns1:n1m) -
            log(ordered1)[ns1p:n1] )
    bid1 = c( cumsum((log(ordered1))^2)[ns1:n1m]/(ns1:n1m) -
            2*cumsum(log(ordered1))[ns1:n1m]*log(ordered1)[ns1p:n1]/(ns1:n1m) +
            ((log(ordered1))^2)[ns1p:n1] )
    dehaan1 = ( 1.0 + bod1 + ( 0.5 / (  bod1^2/bid1 - 1 ) ))
    if (bothTails) {
    n2m = n2-1; ns2 = ns0; ns2p = ns2+1
    bod2 = c( cumsum(log(ordered2))[ns2:n2m]/(ns2:n2m) -
            log(ordered2)[ns2p:n2] )
    bid2 = c(  cumsum((log(ordered2))^2)[ns2:n2m]/(ns2:n2m) -
            2*cumsum(log(ordered2))[ns2:n2m]*log(ordered2)[ns2p:n2]/(ns2:n2m) +
            ((log(ordered2))^2)[ns2p:n2] )
    dehaan2 = ( 1.0 + bod2 + ( 0.5 / (  bod2^2/bid2 - 1 ) )) }

    # PREPARE PLOT:
    y1 = dehaan1[dehaan1 > yrange[1] & dehaan1 < yrange[2]]
    x1 = log10(1:length(dehaan1))[dehaan1 > yrange[1] & dehaan1 < yrange[2]]
    if (bothTails) {
        y2 = dehaan2[dehaan2 > yrange[1] & dehaan2 < yrange[2]]
        x2 = log10(1:length(dehaan2))[dehaan2 > yrange[1] &
            dehaan2 < yrange[2]]
    }

    # Labels:
    if (labels) {
        main = "Deckers - Einmahl - de Haan Estimator"
        xlab = "log scale"
        ylab = "xi"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:
    if (doplot) {
        par(err = -1)
        plot (x1, y1, xlab = xlab, ylab = ylab, ylim = yrange,
            main = main, type = "n")
        if (labels) title(sub = paste("tail depth:", as.character(tail)))
        lines(x1, y1, type = "p", pch = 2, col = "steelblue")
        if (bothTails) lines(x2, y2, type = "p", pch = 6, col = "brown")
        if (labels) grid()
    }

    # CALCULATE INVERSE XI:
    my1 = mean(y1, na.rm = TRUE)
    if (bothTails) my2 = mean(y2, na.rm = TRUE)
    sy1 = sqrt(var(y1, na.rm = TRUE))
    if (bothTails) sy2 = sqrt(var(y2, na.rm = TRUE))
    if (doplot) {
        par(err = -1)
        lines(c(x1[1], x1[length(x1)]), c(my1,my1), type = "l",
            lty = 1, col = "steelblue")
        if (bothTails) lines(c(x2[1], x2[length(x2)]), c(my2, my2),
            type = "l", lty = 1, col = "brown")
    }

    # Result:
    result = list(xi = c(my1, sy1))
    if (bothTails) result = list(xi = c(my1, sy1, my2, sy2))

    # Return Result:
    result
}


################################################################################
