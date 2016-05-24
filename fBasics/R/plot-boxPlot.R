
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
# You should have received A copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA


################################################################################
# FUNCTION:                DESCRIPTION:
#  boxPlot                  Produces a side-by-side standard box plot
#  boxPercentilePlot        Produces a side-by-side box-percentile plot
################################################################################


boxPlot <-
function(x, col = "steelblue", title = TRUE, ...) 
{   
    # A function Implemented by Diethelm Wuertz

    # Description:
    #   Produces a standard box plot
    
    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries' 
    #       or any other object which can be transformed by the function
    #       'as.matrix()' into an object of class 'matrix'.
    
    # Optional Arguments:
    #   las, oma - allows to change style of X labels and creates 
    #       required space below plot. 
    #       Try: e.g. las = 3, and oma = c(9, 0, 0, 0)
    
    
    # FUNCTION:
    
    # Settings:
    x = as.matrix(x)
    assetNames = colnames(x)
    
    # Plot:
    ans = boxplot(as.data.frame(x), col = col, ...)
    abline(h = 0 , lty = 3)
    
    # Add Title:
    if (title) {
        title(main = "Box Plot", ylab = "Value")
    }
    
    # Result:
    colnames(ans$stats) = ans$names
    rownames(ans$stats) = c("lower whisker", "lower hinge", "median", 
        "upper hinge", "upper whisker")
    
    # Return Value:
    invisible(ans)
}   


# ------------------------------------------------------------------------------


boxPercentilePlot <-  
function(x, col = "steelblue", title = TRUE, ...) 
{   
    # A modified copy from Hmisc

    # Description:
    #   Produces a side-by-side box-percentile plot
    
    # Details:
    #   Box-percentile plots are similiar to boxplots, except box-percentile 
    #   plots supply more information about the univariate distributions. At 
    #   any height the width of the irregular "box" is proportional to the 
    #   percentile of that height, up to the 50th percentile, and above the 
    #   50th percentile the width is proportional to 100 minus the percentile. 
    #   Thus, the width at any given height is proportional to the percent of 
    #   observations that are more extreme in that direction. As in boxplots, 
    #   the median, 25th and 75th percentiles are marked with line segments 
    #   across the box. [Source: Hmisc]
    
    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries' 
    #       or any other object which can be transformed by the function
    #       'as.matrix()' into an object of class 'matrix'.
    
    # FUNCTION:
    
    # Settings:
    x = as.matrix(x)
    assetNames = colnames(x)
    n = ncol(x)
    all.x = list()
    for (i in 1:n) all.x[[i]] = as.vector(x[, i])
    centers = seq(from = 0, by = 1.2, length = n)
    ymax = max(sapply(all.x, max, na.rm = TRUE))
    ymin = min(sapply(all.x, min, na.rm = TRUE))
    xmax = max(centers) + 0.5
    xmin = -0.5
    
    # Plot:
    if (length(col) == 1) col = rep(col, times = n)
    plot(c(xmin, xmax), c(ymin, ymax), type = "n",  
        xlab = "", ylab = "", xaxt = "n", ...)
    xpos = NULL
    for (i in 1:n) {
        # plot.values = .bpxAssetsPlot(all.x[[i]], centers[i])
        y = all.x[[i]]
        offset = centers[i]
        y = y[!is.na(y)]
        n = length(y)
        delta = 1/(n + 1)
        prob = seq(delta, 1 - delta, delta)
        quan = sort(y)
        med = median(y)
        q1 = median(y[y < med])
        q3 = median(y[y > med])
        first.half.p = prob[quan <= med]
        second.half.p = 1 - prob[quan > med]
        plotx = c(first.half.p, second.half.p)
        options(warn = -1)
        qx = approx(quan, plotx, xout = q1)$y
        q1.x = c(-qx, qx) + offset
        qx = approx(quan, plotx, xout = q3)$y
        options(warn = 0)
        q3.x = c(-qx, qx) + offset
        q1.y = c(q1, q1)
        q3.y = c(q3, q3)
        med.x = c(-max(first.half.p), max(first.half.p)) + offset
        med.y = c(med, med)
        plot.values = list(x1 = (-plotx) + offset, y1 = quan, x2 = plotx + 
            offset, y2 = quan, q1.y = q1.y, q1.x = q1.x, q3.y = q3.y, 
            q3.x = q3.x, med.y = med.y, med.x = med.x)
        # Continue:
        xpos = c(xpos, mean(plot.values$med.x))
        x.p = c(plot.values$x1, plot.values$x2)
        y.p = c(plot.values$y1, plot.values$y2)
        polygon(x.p, y.p, col = col[i], border = "grey", ...)
        lines(plot.values$x1, plot.values$y1)
        lines(plot.values$x2, plot.values$y2)
        lines(plot.values$q1.x, plot.values$q1.y)
        lines(plot.values$q3.x, plot.values$q3.y)
        lines(plot.values$med.x, plot.values$med.y) 
    }
    axis(side = 1, at = xpos, labels = assetNames, ...)
    abline(h = 0, lty = 3, col = "black")
    
    # Add Title:
    if (title) {
        title(main = "Box Percentiles", ylab = "Value")
    }
   
    # Return Value:
    invisible()
}


################################################################################

