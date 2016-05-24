
# This library is free software, you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation, either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library, if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                   DESCRIPTION:
#  covEllipsesPlot             Displays a covariance ellipses plot
################################################################################


covEllipsesPlot <- 
    function(x = list(), ...)
{
    # Description:
    #   Displays a covariance ellipses plot
    
    # Arguments:
    #   x = a list of at least two covariance matrices   
    
    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   Cov = cov(x); robustCov = assetsMeanCov(x, "MCD")$Sigma
    #   covEllipsesPlot(list(Cov, robustCov))  

    # Source:
    #   Partly based on function covfmEllipsesPlot() from
    #   Package: robust 0.2-2, 2006-03-24
    #   Maintainer: Kjell Konis <konis@stats.ox.ac.uk>
    #   Description: A package of robust methods.
    #   License: Insightful Robust Library License (see license.txt)
    
    # FUNCTION:
    
    # Settings:
    if (length(x) == 0) 
        stop("Input must be a list of at least 2 covariance matrices!")
    nModels = length(x)
    p <- dim(x[[1]])[1]

    # Graphics Frame:
    plot(0, 0, xlim = c(0, p+1), ylim = c(0, p+1), type = "n",
         axes = FALSE, xlab = "", ylab = "", ...)
    box()

    # Correlation Ellipses:
    for(k in 1:nModels) {
        s = sqrt(diag(x[[k]]))
        X = x[[k]] / (s %o% s)
        xCenters = matrix(rep(1:p, p), byrow = TRUE, ncol = p)
        yCenters = matrix(rep(p:1, p), ncol = p)
        points = rep((c(0:180, NA) * pi)/90, (p^2 - p) / 2)
        cors = as.vector(rbind(matrix(X[row(X) < col(X)], nrow = 181, 
            ncol = (p^2 - p)/2, byrow = TRUE), rep(NA, (p^2 - p)/2)))
        xs = 0.475 * cos(points + acos(cors)/2) +
            rep(xCenters[row(xCenters) < col(xCenters)], each = 182)
        ys = 0.475 * cos(points - acos(cors)/2) +
            rep(yCenters[row(xCenters) < col(xCenters)], each = 182)   
        polygon(x = xs, y = ys, density = 0, col = k)
        shift = max(0.2, (p - 8)/88 + 0.2)
        xs = xCenters[row(xCenters) > col(xCenters)]
        ys = yCenters[row(yCenters) > col(yCenters)]
        cors = X[row(X) > col(X)]
        text(xs, ys + (((shift*(nModels - 1))/2) - shift*(k - 1)),
            labels = round(cors, digits = max(1, floor(20/p))),
            col = k, cex = min(1, 90/(p^2)))
    }

    # Diagonal Line:
    lines(c(0.5, p+0.5), c(p+0.5, 0.5), lwd = 2)

    # Correlation - Text:
    text(x = cbind(1:p, rep(p + 0.7, p)), 
        labels = dimnames(X)[[2]], cex = 1, adj = 0)
    text(x = cbind(rep(0.5, p), p:1), 
        labels = dimnames(X)[[1]], cex = 1, adj = 1)
    legend(x = (p+1)/2, y = 0, legend = unlist(paste("-", names(x), "-")), 
        xjust = 0.5, yjust = 0, text.col = 1:nModels, bty = "n")

    # Return Value:
    invisible()
}


################################################################################

