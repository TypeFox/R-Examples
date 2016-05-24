
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
# FUNCTION:             DESCRIPTION:
#  .cov.arw              Energy test for multivariate normality
################################################################################


# Rmetrics:
#   Note that mvoutlier is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 


# Package: mvoutlier                                                     
# Version: 1.4                                                           
# Date: 2009-01-21                                                       
# Title: Multivariate outlier detection based on robust methods          
# Author: Moritz Gschwandtner <e0125439@student.tuwien.ac.at> and        
# 	Peter Filzmoser <P.Filzmoser@tuwien.ac.at>                            
# Maintainer: Peter Filzmoser <P.Filzmoser@tuwien.ac.at>                 
# Depends: R (>= 1.9.0), robustbase, stats                               
# Description: This packages was made for multivariate outlier detection.
# License: GPL version 2 or newer                                        
# URL: http://www.statistik.tuwien.ac.at/public/filz/                    


# ------------------------------------------------------------------------------


.cov.arw <-  
function(x, center, cov, alpha = 0.025, pcrit = NULL)
{
    # Description:
    #   Adaptive reweighted estimator for multivariate location and 
    #   scatter with hard-rejection weights and delta = chi2inv(1-d, p)

    # Arguments
    #   x - Dataset (n x p)
    #   center - Initial location estimator (1 x p)
    #   cov - Initial scatter estimator (p x p)
    #   alpha - Maximum thresholding proportion
    #           (optional scalar, default: alpha = 0.025)
    #   pcrit - critical value for outlier probability
    #           (optional scalar, default values from simulations)

    # Value:
    #   center - Adaptive location estimator (p x 1)
    #   cov - Adaptive scatter estimator (p x p)
    #   cn - Adaptive threshold (scalar)
    #   w - Weight vector (n x 1)
    
    # FUNCTION:
      
    # Settings: 
    x <- getDataPart(x)
    n <- nrow(x)
    p <- ncol(x)
    
    # Critical value for outlier probability based on 
    #   simulations for alpha = 0.025
    if (missing(pcrit)) {
        if (p <= 10) pcrit <- (0.24 - 0.003 * p)/sqrt(n)
        if (p > 10) pcrit <- (0.252 - 0.0018 * p)/sqrt(n)
    }
       
    # Critical value for outlier probability based on 
    #   simulations for alpha = 0.025
    if (p <= 10) pcrit <- (0.24-0.003*p)/sqrt(n)
    if (p > 10) pcrit <- (0.252-0.0018*p)/sqrt(n)
    delta <- qchisq(1 - alpha, p)
    
    d2 <- mahalanobis(x, center, cov)
    d2ord <- sort(d2)
    dif <- pchisq(d2ord,p) - (0.5:n)/n
    i <- (d2ord >= delta) & (dif > 0)
    
    if (sum(i) == 0) 
        alfan <- 0 else alfan <- max(dif[i])
    if (alfan < pcrit) 
        alfan <- 0
    if (alfan > 0) 
        cn <- max(d2ord[n-ceiling(n*alfan)], delta) 
    else 
        cn <- Inf
    w <- d2 < cn
    
    if(sum(w) != 0) {
        center <- apply(x[w, ], 2, mean)
        c1 <- as.matrix(x - rep(1, n) %*% t(center))
        cov <- (t(c1) %*% diag(w) %*% c1) / sum(w)
    }
    
    # Return Value:
    list(center = center, cov = cov, cn = cn, w = w)
}


################################################################################

