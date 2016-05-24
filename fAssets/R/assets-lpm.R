
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
#  assetsLPM                   Computes asymmetric lower partial moments
#  assetsSLPM                  Computes symmetric lower partial moments
################################################################################


assetsLPM <-  
    function(x, tau=colMeans(x), a=1.5, ...)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes LPM and CLPM from multivariate time series
    
    # Arguments:
    #   x - a multivariate time series, a data frame, or any other
    #       rectangular object of assets which can be converted into
    #       a matrix by the function 'as.matrix'. Optional Dates are 
    #       rownames, instrument names are column names.
    
    # References:
    #   Nawrocki 1991, Optimal Algorithms and Lower Partial Moment:
    #       Ex-Post Results
    #   Lee 2006, The Strengths and Limitations of Risk Measures in 
    #       Real Estate: A Review  
    
    # Note:
    #   The output of this function can be used for portfolio
    #   optimization. LPM stands for lower partial moments.
    
    # Example:
    #   LPP <- as.timeSeries(data(LPP2005REC))[, 1:6]; assetsLPM(LPP)
    
    # FUNCTION:
    
    # Transform Input:
    x.mat <- as.matrix(x)
    nCol <- ncol(x)
    nRow <- nrow(x)
    Tau <- matrix(rep(tau, nRow), byrow = TRUE, ncol = nCol)
    TauX <- Tau-x
    X.max <- ((TauX) + abs(TauX))/2
    
    # Compute Lower Partial Moments:
    LPM <- colMeans(X.max^a)
    
    # Compute co-LPMs:
    CLPM <- diag(0, nCol)
    if (a > 1) {
        for (i in 1:nCol) {
            for (j in 1:nCol) {
                CLPM[i, j] <- mean( (X.max[, i])^(a-1) * TauX[, j] )
            }
            CLPM[i, i] <- LPM[i]
        }
    } else if (a == 1) {
        for (i in 1:nCol) {
            for (j in 1:nCol) {
                CLPM[i, j] <- mean( sign( X.max[, i]) * TauX[, j] )
            }
            CLPM[i, i] <- LPM[i]
        }
    }
    
    # Result:
    ans <- list(mu = LPM, Sigma = CLPM)
    attr(ans, "control") <- c(a = a, tau = tau)
        
    # Return Value:
    ans
}


################################################################################
 

assetsSLPM <-  
    function(x, tau=colMeans(x), a=1.5, ...)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes LPM and SLPM from multivariate time series
    
    # Arguments:
    #   x - a multivariate time series, a data frame, or any other
    #       rectangular object of assets which can be converted into
    #       a matrix by the function 'as.matrix'. Optional Dates are 
    #       rownames, instrument names are column names.
    
    # References:
    #   Nawrocki 1991, Optimal Algorithms and Lower Partial Moment:
    #       Ex-Post Results
    #   Lee 2006, The Strengths and Limitations of Risk Measures in 
    #       Real Estate: A Review     
    
    # Note:
    #   The output of this function can be used for portfolio
    #   optimization. SLPM stands for symmetric lower partial moments.
    
    # Example:
    #   LPP = as.timeSeries(data(LPP2005REC))[, 1:6]; assetsSLPM(LPP)
    
    # FUNCTION:
    
    # Transform Input:
    x.mat <- as.matrix(x)
    nCol <- ncol(x)
    nRow <- nrow(x)
    Tau <- matrix(rep(tau, nRow), byrow = TRUE, ncol = nCol)
    TauX <- Tau-x
    X.max <- ((TauX) + abs(TauX))/2
    
    # Compute Lower Partial Moments:
    LPM <- colMeans(X.max^a)
    
    # Compute co-SLPMs:
    SLPM <- LPM^(1/a) %o% LPM^(1/a) * cor(x.mat)
    
    # Result:
    ans <- list(mu = LPM, Sigma = SLPM)
    attr(ans, "control") <- c(a = a, tau = tau)
        
    # Return Value:
    ans
}


################################################################################

    