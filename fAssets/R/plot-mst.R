
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
#  assetsTreePlot              Displays a minimum spanning tree of assets
################################################################################


assetsTreePlot <-
    function(x, labels = TRUE, title = TRUE, box = TRUE,
    method = "euclidian", seed = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a minimum spanning tree of assets
    
    # Arguments:
    #   x -
    #   labels -
    #   title -
    #   box - 
    #   method -
    #   seed -
    
    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   .assetsTreePlot(x) # try several radom choices
    #   .assetsTreePlot(x)
    #   .assetsTreePlot(x) 

    # FUNCTION:

    # Settings:
    if (title) {
        Main = substitute(x)
    } else {
        Main = ""
    }

    # Compute Distance Matrix:
    Order = NULL
    if (class(x) == "dist") {
        DIST = x
    } else {
        # Rank Seed:
        x = series(x)
        if (is.null(seed)) {
            Order = sample(1:ncol(x))
            x = x[, Order]
        }
        DIST = dist(t(x), method[1])
    }
    method = attr(DIST, "method")

    # Compute Minimum Spanning Tree"
    MST = .mst(DIST)

    # Plot Tree:
    .mstPlot(MST, ".nsca", main = Main, ...)
    mtext(paste("Distance Method:", method),
        side = 4, line = 0.1, adj = 0, col = "darkgrey", cex = 0.7)

    # Return Value:
    invisible(list(mst = MST, dist = DIST, order = Order))
}


################################################################################

