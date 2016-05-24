
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
#  assetsFit                   Fits the parameters of a set of assets
#  assetsSim                   Simulates a set of assets
################################################################################


assetsFit <-
function(x, method = c("st", "sn", "sc"), 
  title = NULL, description = NULL, fixed.df = NA, ...)
{
  # A function implemented by Diethelm Wuertz
  
  # Description:
  #   Fits the parameters of a multivariate data set of assets
  #   and returns a list with the values for the mean, the covariance,
  #   the skewness, and the fatness of tails.
  
  # Arguments:
  #   x - A multivariate time series, a data frame, or any other
  #       rectangular object of assets which can be converted into
  #       a matrix by the function as.matrix. Optional Dates are
  #       rownames, instrument names are column names.
  #   method - Which type of distribution should be fitted?
  #       a) st - multivariate skew Student-t
  #       b) sn - multivariate skew Normal
  #       b) sc - multivariate skew-Cauchy
  #       
  
  # Value:
  #   The function returns a list with the following entries:
  #   mu - Mean values of each asset time series
  #   Omega - Covariance matrix of assets
  #   alpha - Skewness vector
  #   df - Degrees of freedom, measures kurtosis
  
  # Notes:
  #   Requires function "msn.mle" and "mst.mle" from R's GPL licensed
  #     contributed package "sn", (C) 1998-2004 A. Azzalini.
  #   The m[method]Fit functions where the "sn" functionality is used
  #     are implemented within the fMultivar package.
  #   The object returned by this function can serve as input for the
  #     function assetsSim().
  
  # FUNCTION:
  
  # Settings:
  assets <- as.matrix(x)
  method <- match.arg(method)
  colNames <- colnames(x)
  
  # Select Distribution:
  FUN <- get(paste0("m", method, "Fit"))
  
  # Fit Parameters:
  ans <- FUN(x=x, trace=FALSE, ...)
  
  # Return Value:
  ans
}

################################################################################

assetsSim <- 
  function(n, method=c("st", "sn", "sc"),
    model=list(beta=rep(0, 2), Omega=diag(2), alpha=rep(0, 2), nu=4), 
    assetNames=NULL) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION
    
    # Match Method:
    method <- match.arg(method)
    
    # Extract Parameters:
    xi <- as.vector(model$beta)
    Omega <- model$Omega
    alpha <- model$alpha
    nu <- model$nu
    
    # Create Random Numbers:
    if (method == "st") ans <- sn::rmst(n, xi, Omega, alpha, nu)
    if (method == "sn") ans <- sn::rmsn(n, xi, Omega, alpha)
    if (method == "sc") ans <- sn::rmsc(n, xi, Omega, alpha)
    
    # Add Optional Asset Names:
    colnames(ans) <- assetNames
    
    # Return Value:
    ans
}


###############################################################################


