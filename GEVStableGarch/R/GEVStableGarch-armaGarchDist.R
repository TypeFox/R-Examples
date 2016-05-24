
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

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
#  FUNCTION:               DESCRIPTION:
#
#  .armaGarchDist          Computes density values for several 
#                          conditional distributions
################################################################################




# ------------------------------------------------------------------------------
# The default is to use the dstable density from package stabledist
# If the user has the 'stable' package from JP Nolan, then the density
# is computed with the dstable.quick function.
# The configuration of the .GSgarch.dstable function was done 
# inside the .onAttach function, but we get some warnings concerning
# .armaGarchDist: no visible global function definition for '.GSgarch.dstable'.
# Therefore, we decided to define the function here and to decide wether 
# or not to use package 'stable' only when function .armaGarchDist is called.
# This is bad for efficiency and we must find another reasonable way to define it


.GSgarch.dstable <- function(x,alpha = 1.5, beta = 0, gamma = 1, 
                             delta = 0, param = 1)
{
  return(stabledist::dstable(x, alpha, beta, gamma, 
                             delta, pm = param))
}



# ------------------------------------------------------------------------------



.armaGarchDist <- 
    function(z, hh, shape = 1.5, skew = 0, 
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
    TOLG = 1e-8) 
{
    # Description:
    #   Calculates the likelihood function for a vector of points (z)
    #   according to the specified distribution.
    #   REMARKS: We choose not to avoid calling this function with out of bound parameters
    #   because in these cases we simply return a huge likelihood to guide 
    #   the optimization algorithm.
    
    # Arguments:
    #   z - vector of points to calculate the llh.
    #   h - vector of conditional variances ????? BETTER DESCRIPTION NEEDED.
    #   cond.dist - name of the conditional distribution
    #   TOLG - general tolerance for arma-garch parameters. 
    #   In the beggining it was set to 1e-5
      
    # FUNCTION:  
      
    # Error treatment of input parameters
    cond.dist = match.arg(cond.dist)
    if(length(z) != length(hh))
        stop("Error: Vectors 'z' and 'hh' have different length.")
    if(sum(is.na(hh)) > 0 || min(hh) == 0)
    {
        return(1e99)
    }
    
    # normal conditional distribution
    if(cond.dist == "norm")
        return(-sum(log(dnorm(x = z/hh)/hh)))
    
    # t-student conditional distribution.
    if(cond.dist == "std")
    {
        if(!(shape > 2))
            return(Inf)
        return(-sum(log(dstd(x = z/hh, nu = shape)/hh)))
    }
    
    # skew t-student conditional (standardized version defined in Wurtz)
    if(cond.dist == "sstd")
    {
        if(!(shape > 2) || !(skew > 0))
            return(1e99)       
        return(-sum(log(dsstd(x = z/hh, nu = shape, xi = skew)/hh)))        

    }
    
    # skew t-student from Fernandez, C. and Steel, M. F. J. (1998)
    if(cond.dist == "skstd")
    {
      if(!(shape > 2) || !(skew > 0))
          return(1e99)     
      return(-sum(log(dskstd(x = z/hh, nu = shape, xi = skew)/hh)))        
    
    }
    
    # GAt distribution
    if(cond.dist == "gat")
    {
        if(!(shape[1] > 0) || !(shape[2] > 0) || !(skew > 0))
           return(1e99)    
        return(-sum(log(dgat(x = z/hh, nu = shape[1], d = shape[2], xi = skew)/hh)))        
    }
    
    # GED conditional distribution.
    if(cond.dist == "ged")
    {
      if(!(shape > 0))
          return(1e99)
      return(-sum(log(dged(x = z/hh, nu = shape)/hh)))
    }
    
    # GEV conditional distribution
    if(cond.dist == "gev")
    {
        # to ensure good mle properties and finiteness of the variance
        # we require shape > -0.5 ( See Jondeau et al. (XXX))
        if( (shape[1] <= -0.5) || (shape[1] >= 0.5)) 
            return(1e99)
        
        # compute density normally with shape != from zero. There was some numerical problems
        # while using function dgev from package fExtremes and thus, we decided to compute the
        # density directly, asumming that the shape parameter is different from zero. Aditionally,
        # we cannot enforce the computation with the shape parameter equal to zero, since the
        # algorithm is not able to move for the optimum in most cases.
        
        # when the shape parameter is small
        if( abs(shape[1]) < 1e-7 )
        {
            llh <- sum(log(hh)) + sum(z/hh) + sum(exp(-z/hh))
            return(llh)
        }
        
        # when the shape parameter is next to zero
        sig <- hh
        y <- 1 + shape * (z/hh)
        llh <- sum(log(sig)) + sum(y^(-1/shape)) + sum(log(y))*(1/shape + 1)
            return(llh)
    }
    
    # stable conditional distribution
    if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
    {
        # Compute density with package 'stable' if it is available
        if(getOption('.stableIsLoaded', default = FALSE) == TRUE)
        {
          .GSgarch.dstable <- function(x,alpha = 1.5, beta = 0, gamma = 1, 
                                       delta = 0, param = 1)
          {
            return(stable::dstable.quick(x, alpha, beta, gamma, 
                                         delta, param))
          }
        }
      
      
      
      
        # Return Big Values if we are out of parameter space

        if( !(shape > 1) || !(shape < 2) || !(abs(skew) < 1))
            return(1e99)        

        y <- z/hh
        if(sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y)))
            return(1e99)
        
        # Compute density either using "stable" or "stabledist" package
        
        if(cond.dist == "stableS0")
            result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape,
                      beta = skew, param = 0)/hh))
        
        if(cond.dist == "stableS1")
          result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape,
                                        beta = skew, param = 1)/hh))
        
        if(cond.dist == "stableS2")
          result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape,
                                        beta = skew, param = 2)/hh))
        
        # Result 
        return(result)
    }
}






################################################################################
