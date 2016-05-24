
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


################################################################################
# FUNCTION:             DESCRIPTION:
#   .col.corrgram
#   .panel.pts
#   .panel.pie
#   .panel.shade
#   .panel.txt
#   .panel.ellipse
#   .corrgram
################################################################################


# Rmetrics:
#   Note that corrgram is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 


# Package: corrgram
# Type: Package
# Title: Plot a correlogram
# Version: 0.1
# Date: 2006-11-28
# Author: Kevin Wright
# Maintainer: Kevin Wright, <kw.statr@gmail.com>
# Description: 
#     Calculates correlation of variables and displays the results graphically.
# License: GPL version 2 or later.
# Packaged: Thu Nov 30


# ------------------------------------------------------------------------------


.col.corrgram = 
function(ncol)
{
    # Colors to use for the corrgram
    # Red > White > Blue
    
    # colorRampPalette(c("red","salmon","white","royalblue","navy"))(ncol)
    
    # colorRampPalette(
    #     c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk"))(ncol)
    
    # heat.colors(ncol)
    
    cm.colors(ncol)
}


# ------------------------------------------------------------------------------


.panel.pts = 
function(x, y, ...)
{
    plot.xy(xy.coords(x, y), type = "p", ...)
    
    box(col = "lightgray")
}


# ------------------------------------------------------------------------------


.panel.pie = 
function(x, y, ...)
{
    # box(col="gray70")
    
    # Coordinates of box
    usr = par()$usr
    minx = usr[1] #min(x, na.rm=TRUE)
    maxx = usr[2] #max(x, na.rm=TRUE)
    miny = usr[3] #min(y, na.rm=TRUE)
    maxy = usr[4] #max(y, na.rm=TRUE)
    
    # Multiply the radius by .97 so the circles do not overlap
    rx = (maxx-minx)/2 * .97
    ry = (maxy-miny)/2 * .97
    centerx = (minx+maxx)/2
    centery = (miny+maxy)/2
    
    segments = 60
    angles = seq(0,2*pi,length=segments)
    circ = cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
    lines(circ[,1], circ[,2], col = 'gray30',...)
    
    # Overlay a colored polygon
    corr = cor(x, y, use = 'pair')
    ncol = 14
    pal = .col.corrgram(ncol)
    col.ind = round(ncol*(corr+1)/2)
    col.pie = pal[col.ind]
    
    segments = round(60*abs(corr),0) # Watch out for the case with 0 segments
    if(segments > 0){
        angles = seq(pi/2, pi/2+(2*pi* -corr), length = segments)
        circ = cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
        circ = rbind(circ, c(centerx, centery), circ[1, ])
        polygon(circ[, 1], circ[, 2], col = col.pie)
    }
}


# ------------------------------------------------------------------------------


.panel.shade = 
function(x, y, ...)
{
    r = cor(x, y, use = 'pair')
    ncol = 14
    pal = .col.corrgram(ncol)
    col.ind = round(ncol*(r+1)/2)
    usr = par("usr")
    
    # Solid fill:
    rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], border = NA)
    
    # Add diagonal lines:
    rect(usr[1], usr[3], usr[2], usr[4], density = 5,
         angle = ifelse(r>0, 45, 135), col = "white")
         
    # Boounding box needs to plot on top of the shading, so do it last.
    box(col = 'lightgray')
}


# ------------------------------------------------------------------------------


.panel.hist = 
function(x, y, ...)
{   # A function implemented by Diethelm Wuertz

    # Settings:
    r = cor(x, y, use = 'pair')
    ncol = 14
    pal = .col.corrgram(ncol)
    col.ind = round(ncol*(r+1)/2)
    usr = par("usr")
    
    # Hexagonal Binning:
    object = hexBinning(x, y, bins = 10) 
    X = object$x
    Y = object$y
    rx = min(diff(unique(sort(X))))
    ry = min(diff(unique(sort(Y))))
    rt = 2 * ry
    u = c(rx, 0, -rx, -rx, 0, rx)
    v = c(ry, rt, ry, -ry, -rt, -ry)/3
    N = length(col)
    z = object$z
    zMin = min(z)
    zMax = max(z)
    Z = (z - zMin)/(zMax - zMin)
    Z = trunc(Z * (N - 1) + 1)
    for (i in 1:length(X)) {
        polygon(u + X[i], v + Y[i], col = col[Z[i]], border = "white")
    }
    # points(object$xcm, object$ycm, pch = 19, cex = 1/3, col = "black")
    box(col = 'lightgray')
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.panel.txt = 
function(x = 0.5, y = 0.5, txt, cex, font)
{
    text(x, y, txt, cex = cex, font = font)
    
    # box(col = "lightgray")
}


# ------------------------------------------------------------------------------


.panel.minmax = 
function(x, ...)
{
    # Put the minimum in the lower-left corner and the
    # maximum in the upper-right corner
    minx = round(min(x, na.rm = TRUE),2)
    maxx = round(max(x, na.rm = TRUE),2)
    text(minx, minx, minx, cex = 1, adj = c(0, 0))
    text(maxx, maxx, maxx, cex = 1, adj = c(1, 1))
}


# ------------------------------------------------------------------------------


.panel.ellipse = 
function(x, y, ...)
{
    # Draw an Ellipse:
    dfn = 2
    dfd = length(x)-1
    shape = var(cbind(x,y), na.rm = TRUE)
    keep = (!is.na(x) & !is.na(y))
    center = c(mean(x[keep]),mean(y[keep]))
    radius = sqrt(dfn*qf(.68,dfn,dfd))
    segments = 75
    angles = seq(0,2*pi,length=segments)
    unit.circle = cbind(cos(angles),sin(angles))
    ellipse.pts = t(center+radius*t(unit.circle%*%chol(shape)))
    ellx = ellipse.pts[, 1]
    elly = ellipse.pts[, 2]
    
    # Truncate Ellipse at min/max or at Bounding Box
    usr = par()$usr
    minx = usr[1] #min(x, na.rm=TRUE)
    maxx = usr[2] #max(x, na.rm=TRUE)
    miny = usr[3] #min(y, na.rm=TRUE)
    maxy = usr[4] #max(y, na.rm=TRUE)
    ellx = ifelse(ellx < minx, minx, ellx)
    ellx = ifelse(ellx > maxx, maxx, ellx)
    elly = ifelse(elly < miny, miny, elly)
    elly = ifelse(elly > maxy, maxy, elly)
    # lines(ellx, elly, col = 'gray30', ...)
    
    # Polygon:
    r = cor(x, y, use = 'pair')
    ncol = 14
    pal = .col.corrgram(ncol)
    col.ind = round(ncol*(r+1)/2)
    polygon(ellx, elly, col = pal[col.ind])
    
    # Add a lowess line through the ellipse:
    ok = is.finite(x) & is.finite(y)
    if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = 2/3, iter = 3), col = "red", ...) 
    
    box(col = 'lightgray') 
}


# ------------------------------------------------------------------------------


.panel.copula = 
function (x, y, ...) 
{   # A function Implemented by Diethelm Wuertz
 
    R1 = as.vector(x)
    R2 = as.vector(y)
    
    r1 = R1[R1 != 0 & R2 != 0]
    fit1 = nigFit(r1, doplot = FALSE)
    estim1 = fit1@fit$estimate
    p1 = pnig(r1, estim1[1], estim1[2], estim1[3], estim1[4]) 
    
    r2 = R2[R1 != 0 & R2 != 0]
    fit2 = nigFit(r2, doplot = FALSE)
    estim2 = fit2@fit$estimate
    p2 = pnig(r2, estim2[1], estim2[2], estim2[3], estim2[4]) 
    
    # Rescale to get plotted
    x = (max(r1)-min(r1))*p1 + min(r1)
    y = (max(r2)-min(r2))*p2 + min(r2)

    plot.xy(xy.coords(x, y), type = "p", pch = 19, cex = 0.5, ...)

    box(col = "lightgray")
}


# ----------------------------------------------------------------------------


.corrgram = 
function (x, labels, panel = .panel.shade, ...,
    lower.panel = panel, upper.panel = panel,
    diag.panel = NULL, text.panel = textPanel,
    label.pos = 0.5,
    cex.labels = NULL, font.labels = 1,
    row1attop = TRUE, gap = 1)
{
    textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
            text(x, y, txt, cex = cex, font = font)

    localAxis = function(side, x, y, xpd, bg, col = NULL, main, oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        if(side %% 2 == 1) Axis(x, side = side, xpd = NA, ...)
        else Axis(y, side = side, xpd = NA, ...)
    }

    localPlot = function(..., 
        main, oma, font.main, cex.main) plot(...)
    localLowerPanel = function(..., main, oma, font.main, cex.main)
        lower.panel(...)
    localUpperPanel = function(..., main, oma, font.main, cex.main)
        upper.panel(...)
    localDiagPanel = function(..., main, oma, font.main, cex.main)
        diag.panel(...)

    dots = list(...)
    nmdots = names(dots)
    
    if (!is.matrix(x)) {
        x = as.data.frame(x)
        for(i in seq(along=names(x))) {
            if(is.factor(x[[i]]) || is.logical(x[[i]]))
               x[[i]] = as.numeric(x[[i]])
            if(!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    } else if (!is.numeric(x)) {
        stop("non-numeric argument to 'pairs'")
    }
    
    panel = match.fun(panel)
    if((has.lower = !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel = match.fun(lower.panel)
    if((has.upper = !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel = match.fun(upper.panel)
    if((has.diag  = !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel = match.fun(diag.panel)

    if(row1attop) {
        tmp = lower.panel
        lower.panel = upper.panel
        upper.panel = tmp
        tmp = has.lower
        has.lower = has.upper
        has.upper = tmp
    }

    nc = ncol(x)
    if (nc < 2) stop("only one column in the argument to 'pairs'")
    
    has.labs = TRUE
    if (missing(labels)) {
        labels = colnames(x)
        if (is.null(labels)) labels = paste("var", 1:nc)
    } else if(is.null(labels)) {
        has.labs = FALSE
    }
    oma = if("oma" %in% nmdots) dots$oma else NULL
    main = if("main" %in% nmdots) dots$main else NULL
    if (is.null(oma)) {
        oma = c(4, 4, 4, 4) 
        if (!is.null(main)) oma[3] = 6
    }
    opar = par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))

    for (i in if(row1attop) 1:nc else nc:1)
        for (j in 1:nc) {
            localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ...)
            if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
                mfg = par("mfg")
                if(i == j) {
                    if (has.diag) {
                        localDiagPanel(as.vector(x[, i]), ...)
                    }
                    if (has.labs) {
                        par(usr = c(0, 1, 0, 1))
                        if(is.null(cex.labels)) {
                            l.wid = strwidth(labels, "user")
                            cex.labels = max(0.8, min(2, .9 / max(l.wid)))
                        }
                        text.panel(0.5, label.pos, labels[i],
                                   cex = cex.labels, font = font.labels)
                    }
                } else if(i < j) {
                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                } else {
                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                }
                if (any(par("mfg") != mfg))
                    stop("the 'panel' function made a new plot")
            } else {
                par(new = FALSE)
            }

        }
        if (!is.null(main)) {
            font.main = 
                if("font.main" %in% nmdots) dots$font.main else par("font.main")
            cex.main = 
                if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
            mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}


################################################################################

