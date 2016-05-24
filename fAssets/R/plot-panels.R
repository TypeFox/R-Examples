
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
# FUNCTION:                   DESCRIPTION:
# .txtPanel                    a diagonal text panel
# .minmaxPanel                 a diagonal minmax text panel
# .histPanel                   a diagonal histogram panel
# FUNCTION:                   DESCRIPTION:
# .ptsPanel                    an off-diagonal points panel
# .piePanel                    an off-diagonal pie panel
# .piePtsPanel                 an off-diagonal pie/points panel
# .shadePanel                  an off-diagonal shade panel
# .ellipsePanel                an off-diagonal ellipse panel
# .cortestPanel                an off-diagonal cortest panel
# .lowessPanel                 an off-diagonal lowess panel
# .numberPanel                 an off-diagonal lowess panel
################################################################################


.txtPanel <- 
    function(x = 0.5, y = 0.5, txt, cex, font, col.box = "white")
{  
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    # Text Panel:
    text(x, y, txt, cex = cex, font = font)
    
    # Add Box:
    box(col = col.box)
}


# ------------------------------------------------------------------------------


.minmaxPanel <- 
    function(x, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    # Put the minimum in the lower-left corner and the
    #   maximum in the upper-right corner
    minx <- round(min(x, na.rm = TRUE), 2)
    maxx <- round(max(x, na.rm = TRUE), 2)
    text(minx, minx, minx, cex = 1, adj = c(0,0))
    text(maxx, maxx, maxx, cex = 1, adj = c(1,1))
}


# ------------------------------------------------------------------------------


.histPanel <- 
    function(x, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts 
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


# ------------------------------------------------------------------------------


.ptsPanel <- 
    function(x, y, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    plot.xy(xy.coords(x, y), type = "p", ...)
    box(col = col.box)
}


# ------------------------------------------------------------------------------


.piePanel <- 
    function(x, y, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    # Coordinates of box:
    usr <- par()$usr
    minx <- usr[1] #min(x, na.rm = TRUE)
    maxx <- usr[2] #max(x, na.rm = TRUE)
    miny <- usr[3] #min(y, na.rm = TRUE)
    maxy <- usr[4] #max(y, na.rm = TRUE)
    # Multiply the radius by 0.97 so the circles do not overlap
    rx <- (maxx-minx)/2 * 0.97
    ry <- (maxy-miny)/2 * 0.97
    centerx <- (minx+maxx)/2
    centery <- (miny+maxy)/2
    
    segments <- 60
    angles <- seq(0, 2*pi, length = segments)
    circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
    lines(circ[,1], circ[,2], col = 'gray30',...)
    
    # Overlay a colored polygon
    corr <- cor(x, y, use = 'pair')
    ncol <- 14
    pal <- .col.corrgram(ncol)
    col.ind <- round(ncol*(corr+1)/2)
    col.pie <- pal[col.ind]
    
    # Watch out for the case with 0 segments:
    segments <- round(60*abs(corr),0)
    if(segments > 0) {
        angles <- seq(pi/2, pi/2+(2*pi* -corr), length=segments)
        circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
        circ <- rbind(circ, c(centerx, centery), circ[1, ])
        polygon(circ[,1], circ[,2], col = col.pie)
    }
    
    box(col = col.box)
}


# ------------------------------------------------------------------------------


.piePtsPanel <-
    function(x, y, col.box = "white", ...)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Example:
    #   x = series(100*as.timeSeries(data(LPP2005REC))[, 1:6])
    #   pairs(x, tick = FALSE)
    
    # FUNCTION:
    
    # Pie Panel:
    .piePanel(x, y, col.box = "white", ...)  
    
    # Add Points:
    points(x, y, ...)
}


# ------------------------------------------------------------------------------


.shadePanel <- 
    function(x, y, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    r <- cor(x, y, use = 'pair')
    ncol <- 14
    pal <- .col.corrgram(ncol)
    col.ind <- round(ncol*(r+1)/2)
    usr <- par("usr")
    
    # Solid fill:
    rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], border = NA)
    
    # Add diagonal lines:
    rect(usr[1], usr[3], usr[2], usr[4], density = 5,
         angle = ifelse(r>0, 45, 135), col="white")
         
    # Bounding box needs to plot on top of the shading, so do it last.
    
    box(col = col.box)
}


# ------------------------------------------------------------------------------


.ellipsePanel <- 
    function(x, y, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    # Draw an ellipse:
    #  box(col="white")
    dfn <- 2
    dfd <- length(x)-1
    shape <- var(cbind(x,y),na.rm=TRUE)
    keep <- (!is.na(x) & !is.na(y))
    center <- c(mean(x[keep]),mean(y[keep]))
    radius <- sqrt(dfn*qf(.68,dfn,dfd))
    segments <- 75
    angles <- seq(0,2*pi,length=segments)
    unit.circle <- cbind(cos(angles),sin(angles))
    ellipse.pts <- t(center+radius*t(unit.circle%*%chol(shape)))
    ellx <- ellipse.pts[,1]
    elly <- ellipse.pts[,2]
    # Truncate ellipse at min/max or at bounding box
    usr <- par()$usr
    minx <- usr[1] #min(x, na.rm=TRUE)
    maxx <- usr[2] #max(x, na.rm=TRUE)
    miny <- usr[3] #min(y, na.rm=TRUE)
    maxy <- usr[4] #max(y, na.rm=TRUE)
    ellx <- ifelse(ellx < minx, minx, ellx)
    ellx <- ifelse(ellx > maxx, maxx, ellx)
    elly <- ifelse(elly < miny, miny, elly)
    elly <- ifelse(elly > maxy, maxy, elly)
    lines(ellx, elly, col='gray30',...)
    
    # Fill Ellipse:
    # polygon(ellx, elly, col = "blue", ...)
    
    # Add a lowess line through the ellipse
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = 2/3, iter = 3), col = "red", ...) 
    
    box(col = col.box) 
}


# ------------------------------------------------------------------------------


.cortestPanel <-
    function(x, y, cex, col, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    if (missing(col)) col = NULL
    usr = par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = abs(cor(x, y))
    txt = format(c(r, 0.123456789), digits = 3)[1]
    test = cor.test(x, y)
    Signif = symnum(test$p.value, corr = FALSE, na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("*** ", "** ", "* ", ". ", "  "))
    text(0.5, 0.5, txt, cex = 1, col = NULL, ...)
    text(0.8, 0.8, Signif, cex = 1.5, col = col, ...)
}


# ------------------------------------------------------------------------------


.lowessPanel =
    function (x, y, col.box = "white", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    points(x, y, ...)
    ok = is.finite(x) & is.finite(y)
    if (any(ok)) lines(lowess(x[ok], y[ok]), col = "brown")
    
    box(col = col.box)
}


# ------------------------------------------------------------------------------


.numberPanel <-
    function(x, y, cex, col, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    if (missing(col)) col = NULL
    usr = par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    
    # Correletion Coefficient
    number = as.character(round(100*cor(x, y)))
    
    text(0.5, 0.5, number, cex = 1, col = NULL, ...)
}


################################################################################

