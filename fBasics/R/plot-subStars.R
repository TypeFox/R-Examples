
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
# FUNCTION:                 SLIDER MENU:
#  .subStars                 Creates a sub stars plot
################################################################################


.subStars <- 
function(x, keyOffset = c(0, 0), loc = NULL, palette = NULL) 
{
    # Description:
    #   Creates a sub stars plot
    
    # Details:
    #   this is the stars::graphics() function modified in such a
    #   way, that it works together with the subplot builtin.
    
    # FUNCTION:
    
    # Palette:
    if (is.null(palette)) palette = palette(rainbow(12, s = 0.6, v = 0.75))
    
    # Graph Frame:
    xCol = ncol(x)
    yCol = nrow(x)
    NY = NX = ceiling(sqrt(xCol))
    if (NX * NY == xCol) NY = NY + 1
    if (is.null(loc)) {
        for (nx in 1:NY) for (ny in 1:NX) loc = rbind(loc, c(nx, ny))
        loc = loc[1:xCol, ]
        loc[, 2] = NY + 1 - loc[, 2]
    }  
   
    # Stars - Settings: 
    x = t(x)
    locations = loc 
    len = 0.4
    key.loc = c(NX + 1, 1) + keyOffset
    xlim = c(1, NX + 0.5) 
    ylim = c(0, NY + 1)
        
    op <- par("xpd")
    on.exit(par(op))
    par(xpd = TRUE)
    
    labels = dimnames(x)[[1]]
    key.labels = dimnames(x)[[2]]
    
    n.loc <- nrow(x)
    n.seg <- ncol(x)

    xloc <- locations[, 1]
    yloc <- locations[, 2]
    
    angles <- seq.int(0, 2 * pi, length.out = n.seg + 1)[-(n.seg + 1)]
    x <- apply(x, 2, 
    function(x) (x - min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE)))
    
    mx <- max(x <- x * len)
    if (is.null(xlim)) xlim <- range(xloc) + c(-mx, mx)
    if (is.null(ylim)) ylim <- range(yloc) + c(-mx, mx)
    deg <- pi/180
    
    plot(0, type = "n",  xlim = xlim, ylim = ylim, main = "", 
         sub = "", xlab = "", ylab = "", asp = 1, axes = FALSE)

    s.x <- xloc + x * rep.int(cos(angles), rep.int(n.loc, n.seg))
    s.y <- yloc + x * rep.int(sin(angles), rep.int(n.loc, n.seg))
    
    aangl <- c(angles, if (TRUE) 2 * pi else pi)
    for (i in 1:n.loc) {
        px <- py <- numeric()
        for (j in 1:n.seg) {
            k <- seq.int(from = aangl[j], to = aangl[j + 1], by = 1 * deg)
            px <- c(px, xloc[i], s.x[i, j], x[i, j] * cos(k) + xloc[i], NA)
            py <- c(py, yloc[i], s.y[i, j], x[i, j] * sin(k) + yloc[i], NA)
        }
        polygon(px, py, col = 1:n.seg, lwd = 0.25, lty = 1)
    }
    text(xloc, yloc - mx + 0.04, labels, cex = 0.7, font = 2, adj = c(0.5, 1))

    key.x <- len * cos(angles) + key.loc[1]
    key.y <- len * sin(angles) + key.loc[2]
    px <- py <- numeric()
    for (j in 1:n.seg) {
        k <- seq.int(from = aangl[j], to = aangl[j + 1], by = 1 * deg)
        px <- c(px, key.loc[1], key.x[j], len * cos(k) + key.loc[1], NA)
        py <- c(py, key.loc[2], key.y[j], len * sin(k) + key.loc[2], NA)
    }
    polygon(px, py, col = 1:n.seg, lwd = 0.25, lty = 1)
    
    lab.angl <- angles + (angles[2] - angles[1])/2 
    label.x <- 1.15 * len * cos(lab.angl) + key.loc[1]
    label.y <- 1.15 * len * sin(lab.angl) + key.loc[2]
    
    for (k in 1:n.seg) {
        text.adj <- c(if (lab.angl[k] < 90 * deg || lab.angl[k] > 
            270 * deg) 0 else if (lab.angl[k] > 90 * deg && 
            lab.angl[k] < 270 * deg) 1 else 0.5, if (lab.angl[k] <= 
            90 * deg) (1 - lab.angl[k]/(90 * deg))/2 else if (lab.angl[k] <= 
            270 * deg) (lab.angl[k] - 90 * deg)/(180 * deg) else 1 - 
            (lab.angl[k] - 270 * deg)/(180 * deg))
        text(label.x[k], label.y[k], labels = key.labels[k], 
            cex = 0.7, adj = text.adj)
    }

    # Return Value:
    invisible(locations)
}


################################################################################

