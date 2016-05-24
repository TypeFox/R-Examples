
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
# FUNCTION:                DESCRIPTION:
#  assetsArrange            Rearranges the columns in a deta set of assets
# FUNCTION:                DESCRIPTION:
#  pcaArrange               Returns PCA correlation ordered column names
#  hclustArrange            Returns hierarchical clustered column names 
#  abcArrange               Returns sorted column names 
#  orderArrange             Returns ordered column names 
#  sampleArrange            Returns sampled column names 
#  statsArrange             Returns statistically rearranged column names
################################################################################


assetsArrange <-
    function(x, method = c("pca", "hclust", "abc"), ...)
{ 
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns ordered column names of a time Series

    # Arguments:
    #   x - S4 object of class 'timeSeries'
    
    # FUNCTION:
    
    # Settings:
    method <- match.arg(method)
    FUN <- paste(method, "Arrange", sep = "")
    arrange <- match.fun(FUN)
    
    # Return Value:
    arrange(x, ...)
}
    
# ------------------------------------------------------------------------------


pcaArrange <-
    function(x, robust = FALSE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns PCA correlation ordered column names

    # Arguments:
    #   x - S4 object of class 'timeSeries'
        
    # Notes:
    #   Requires package "robustbase".
    
    # FUNCTION:

    # Order:
    if (robust) {
        x.cor <- robustbase::covMcd(as.matrix(x), cor = TRUE, ...)$cor
    } else {
        x.cor <- cor(as.matrix(x), ...)
    }
    x.eigen <- eigen(x.cor)$vectors[,1:2]
    e1 <- x.eigen[, 1]
    e2 <- x.eigen[, 2]
    Order <- order(ifelse(e1 > 0, atan(e2/e1), atan(e2/e1)+pi))
    ans <- colnames(as.matrix(x))[Order]

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


hclustArrange <- 
    function(x, method = c("euclidean", "complete"), ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns hierarchical clustered column names

    # Arguments:
    #   x - S4 object of class 'timeSeries'
    #   ...
    #       method - the agglomeration method to be used. This should
    #           be (an unambiguous abbreviation of) one of "ward", "single",
    #           "complete", "average", "mcquitty", "median" or "centroid".

    # FUNCTION:

    # Order:
    Order <- hclust(
        dist(t(as.matrix(x)),
        method = method[1]), 
        method = method[2], ...)$order
    ans <- colnames(as.matrix(x))[Order]

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


abcArrange <- 
    function(x, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns sorted column names of a time Series

    # Arguments:
    #   x - S4 object of class 'timeSeries'
    
    # FUNCTION:

    # Sort:
    ans <- sort(colnames(as.matrix(x)), ...)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


orderArrange <-
    function(x, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns ordered column names of a time Series

    # Arguments:
    #   x - S4 object of class 'timeSeries'
    
    # FUNCTION:

    # Order:
    ans <- order(colnames(as.matrix(x)), ...)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


sampleArrange <- 
    function(x, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns sampled column names of a time Series

    # Arguments:
    #   x - S4 object of class 'timeSeries'
    
    # FUNCTION:

    # Sample:
    ans <- sample(colnames(as.matrix(x)), ...)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


statsArrange <- 
    function(x, FUN = colMeans, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns statistically rearranged column names

    # Arguments:
    #   x - S4 object of class 'timeSeries'
    
    # Note:
    #   Example of function Candidates:
    #   colStats        calculates column statistics,
    #   colSums         calculates column sums,
    #   colMeans        calculates column means,
    #   colSds          calculates column standard deviations,
    #   colVars         calculates column variances,
    #   colSkewness     calculates column skewness,
    #   colKurtosis     calculates column kurtosis,
    #   colMaxs         calculates maximum values in each column,
    #   colMins         calculates minimum values in each column,
    #   colProds        computes product of all values in each column,
    #   colQuantiles    computes quantiles of each column.

    # FUNCTION:

    # Apply colStats Function:
    fun <- match.fun(FUN)
    Sort <- sort(fun(x, ...))
    Order <- names(Sort)
    ans <- colnames(as.matrix(x)[, Order])
    attr(ans, "control") <- Sort

    # Return Value:
    ans
}


################################################################################


