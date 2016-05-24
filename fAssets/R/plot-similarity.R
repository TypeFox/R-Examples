
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
#  assetsDendrogramPlot        Displays hierarchical clustering dendrogram
#  assetsCorEigenPlot          Displays ratio of the largest two eigenvalues
################################################################################


assetsDendrogramPlot <-
    function(x, labels = TRUE, title = TRUE, box = TRUE,
    method = c(dist = "euclidian", clust = "complete"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays hierarchical clustering dendrogram

    # FUNCTION:

    # Compute Distance Matrix:
    if (class(x) == "dist") {
        DIST = x
    } else {
        X = t(series(x))
        DIST = dist(X, method[1])
    }

    # Hierarchical Clustering:
    ans = hclust(DIST, method = method[2])

    # Plot Dendrogram:
    if (labels) {
        plot(ans, xlab = "", main = "", sub = "", ...)
        mtext(paste(
            "Distance Method:", method[1], " | ",
            "Clustering Method:", method[2]),
            side = 4, line = 0.1, adj = 0, col = "darkgrey")
    } else {
        plot(ans, ann = FALSE, ...)
    }

    # Add Box:
    if (box) {
        box()
    }

    # Add Optional Title:
    if (title) {
        title(main = "Dendrogram", sub = "", xlab = "", ylab = "Heights")
    }

    # Return Value:
    invisible(list(dist = DIST, hclust = ans))
}


# ------------------------------------------------------------------------------


assetsCorEigenPlot <-
    function(x, labels = TRUE, title = TRUE, box = TRUE,
    method = c("pearson", "kendall", "spearman"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays ratio of the largest two eigenvalues

    # Arguments:
    #   x - a timeSeries object or any other rectangular object
    #       which can be transformed by the function as. matrix
    #       into a numeric matrix.

    # Example:
    #   assetsCorEigenPlot(x=100*as.timeSeries(data(LPP2005REC)))

    # FUNCTION:

    # Settings:
    stopifnot(is.timeSeries(x))
    x = series(x)
    method = match.arg(method)

    # Plot:
    x.cor = cor(x, use = 'pair', method = method)
    x.eig = eigen(x.cor)$vectors[, 1:2]
    e1 = x.eig[, 1]
    e2 = x.eig[, 2]
    plot(e1, e2, col = 'white', ann = FALSE,
        xlim = range(e1, e2), ylim = range(e1, e2), ...)
    abline(h = 0, lty = 3, col = "grey")
    abline(v = 0, lty = 3, col = "grey")
    arrows(0, 0, e1, e2, cex = 0.5, col = "steelblue", length = 0.1)
    text(e1, e2, rownames(x.cor), ...)

    # Labels:
    if (labels) {
        mtext(method, side = 4, adj = 0, cex = 0.7, col = "grey")
    }

    # Add Box:
    if (box) {
        box()
    }

    # Add Title:
    if(title) {
        title(main = "Eigenvalue Ratio Plot", sub = "",
            xlab = "Eigenvalue 1", ylab = "Eigenvalue 2")
    }

    # Return Value:
    invisible()
}


################################################################################

