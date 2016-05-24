
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
#  assetsOutliers              Detects outliers in multivariate assets sets
################################################################################


assetsOutliers <- 
    function (x, center, cov, ...)
{   
    # An adapted copy from contributed R package mvoutlier

    # Description:
    #   Detects outliers in a multivariate set of assets
    
    # Arguments:
    
    # Source:
    #   The code concerned with the outliers is from R package "mvoutliers"
    #   Moritz Gschwandtner <he0125439@student.tuwien.ac.ati>
    #   Peter Filzmoser <hP.Filzmoser@tuwien.ac.ati> 
    
    # References:
    #   P. Filzmoser, R.G. Garrett, and C. Reimann (2005). 
    #   Multivariate Outlier Detection in Exploration Geochemistry. 
    #   Computers & Geosciences.
       
    # FUNCTION:
    
    # Check timeSeries Input:
    stopifnot(is.timeSeries(x))
    tS <- x
    x <- as.matrix(x)
    
    # Critical Values:
    n = nrow(x)
    p = ncol(x)
    if (p <= 10) pcrit = (0.240 - 0.0030 * p)/sqrt(n)
    if (p  > 10) pcrit = (0.252 - 0.0018 * p)/sqrt(n)
    delta <- qchisq(0.975, p)
    
    # Compute Mahalanobis Squared Distances:
    d2 <- mahalanobis(x, center, cov)
    
    # Detect Outliers:
    d2ord <- sort(d2)
    dif <- pchisq(d2ord, p) - (0.5:n)/n
    i <- (d2ord >= delta) & (dif > 0)
    if (sum(i) == 0) alfan = 0 else alfan = max(dif[i])
    if (alfan < pcrit) alfan = 0
    if (alfan > 0) cn = max(d2ord[n-ceiling(n*alfan)], delta) else cn = Inf
    w <- d2 < cn
    m <- apply(x[w, ], 2, mean)
    c1 <- as.matrix(x - rep(1, n) %*% t(m))
    c <- (t(c1) %*% diag(w) %*% c1)/sum(w)
    
    # Identify Outliers:
    outliers <- (1:dim(x)[1])[!w]
    if (length(outliers) == 0) {
        outliers <- NA
    } else {
        names(outliers) <- rownames(x)[outliers]
    }
    
    # Compose Result:
    ans <- list(
        center = m, 
        cov = c, 
        cor = cov2cor(c), 
        quantile = cn, 
        outliers = outliers, 
        series = tS[outliers, ])
    
    # Return Value:
    ans
}


################################################################################

