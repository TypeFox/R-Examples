
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


################################################################################
# FUNCTION:             DESCRIPTION:
#  sampleLmoments        Computes sample L-moments
################################################################################


sampleLmoments <- 
function(x, rmax=4) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description
    #   Computes sample L-moments
    
    # Note:
    #   This function is borrowed from package ...
    #   Author:
    
    # FUNCTION:
    
    # L Moments:
    data <- as.matrix(x)
    n <- dim(data)[1]
    p <- dim(data)[2]
    x <- array(, c(p, n))
    L  <- array(, c(p, rmax))
    for (i in 1:p) x[i, ] <- sort(data[, i])
    if (rmax == 1) return(rowMeans(x))
    bcoef <- array(, c(rmax, n))
    bcoefm <- array(, c(rmax, p, n))
    b <- array(, c(p, rmax))
    bcoef[1, ] <- seq(0, 1, by = (1/(n-1)))
    bcoefm[1, , ] <- t(array(rep(bcoef[1, ], p), c(n, p))) 
    b[, 1] <- rowMeans(x)
    b[, 2] <- rowMeans(bcoefm[1, , ] * x)
    L[, 1] = b[, 1]
    if (rmax > 2) {
        for (r in 2:(rmax-1)) {
            rr <- r+1
            bcoef[r, ]<-bcoef[r-1,]*seq((-(r-1)/(n-r)),1, by = (1/(n-r)))
            bcoefm[r, , ]<-t(array(rep(bcoef[r,],p),c(n, p))) 
            b[, rr] <- rowMeans(bcoefm[r, , ]*x)
        }
    } 
    for (r in 1:(rmax-1)) {
        L[, r+1] <- 0
        for (k in 0:r) {
            kk <- k+1 
            L[, r+1] <- L[, r+1]+(-1)^(r-k)*gamma(r+k+1) / 
                (gamma(k+1)^2) / gamma(r-k+1)*b[, kk] 
        }
    }
    L = as.vector(L)
    names(L) = paste("L", 1:rmax, sep = "")
    
    # Return Value:
    L
}


################################################################################

