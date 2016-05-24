
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
# FUNCTION:                     DESCRIPTION:
#  setType<-                     Sets type of portfolio optimization
#  setOptimize<-                 Sets what to optimze, minRisk or maxRetururn
#  setEstimator<-                Sets name of mean-covariance estimator
#  setTailRisk<-                 Sets tail dependency matrix
#  setParams<-                   Sets optional model parameters
# FUNCTION:                     DESCRIPTION:
#  setWeights<-                  Sets weights vector
#  setTargetReturn<-             Sets target return value
#  setTargetRisk<-               Sets target return value
#  setRiskFreeRate<-             Sets risk-free rate value
#  setNFrontierPoints<-          Sets number of frontier points
#  setStatus<-                   Sets portfolio status information
# FUNCTION:                     DESCRIPTION:
#  setSolver<-                   Sets name of desired solver
#  setObjective<-                Sets objective function name
#  setTrace<-                    Sets solver's trace flag
################################################################################


"setType<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                
    #   Sets the portfolio type for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
    
    # Type ?
    spec@model$type = value
    if (value == "CVaR") setSolver(spec) <- "solveRglpk"
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setOptimize<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                
    #   Sets the portfolio type for a portfolio structure
     
    # Arguments:
    
    # FUNCTION:
    
    # Type ?
    spec@model$optimize = value
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setEstimator<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                  
    #   Sets the type of mean-cov estimator for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
    
    # Estimator ?
    spec@model$estimator = value 
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setParams<-" <- 
    function(spec, name, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                  
    #   Sets optional parameters for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:

    # Extend Parameter List:
    spec@model$params[name] <- value 
    
    # Return Value:
    spec
}


################################################################################


"setWeights<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                    
    #   Sets the weights vector for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
    
    # Weights ?
    spec@portfolio$weights = value
    spec@portfolio$targetReturn = NA
    spec@portfolio$targetRisk = NA
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setTargetReturn<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                                   
    #   Sets the target return value for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:

    # Target Return ?
    spec@portfolio$targetReturn = value
    spec@portfolio$weights = NA
    spec@portfolio$targetRisk = NA
    
    # What to optimize ?
    spec@model$optimize = "minRisk"
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setTargetRisk<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                                   
    #   Sets the target return value for a portfolio structure
  
    # Arguments:
    
    # FUNCTION:
 
    # Target Return ?
    spec@portfolio$targetRisk = value
    spec@portfolio$weights = NA
    spec@portfolio$targetReturn = NA
     
    # What to optimize ?
    spec@model$optimize = "maxReturn"
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setAlpha<-" <-
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                  
    #   Sets the CVaR alpha significance value for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
 
    # Estimator ?
    spec@model$params$alpha = value 
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setRiskFreeRate<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                                   
    #   Sets the risk free rate for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
    
    # Check Validity:
    stopifnot(is.numeric(value))
    stopifnot(length(value) == 1)
    
    # Risk-Free Rate ?
    spec@portfolio$riskFreeRate = value
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setNFrontierPoints<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                                   
    #   Sets the number of frontier points for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
    
    # Check Validity:
    stopifnot(is.numeric(value))
    stopifnot(length(value) == 1)
    stopifnot(value > 0)
    
    # Risk-Free Rate ?
    spec@portfolio$nFrontierPoints = value
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setStatus<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                                   
    #   Sets portfolio status information
    
    # Arguments:
    
    # FUNCTION:
    
    # Check Validity:
    stopifnot(is.numeric(value))
    stopifnot(length(value) == 1)
    
    # Risk-Free Rate ?
    spec@portfolio$status = value
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setTailRisk<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:                    
    #   Sets the tail risk value for a portfolio structure
    
    # Arguments:
    #   value - a list with two matrix elements, $lower and $upper, 
    #       with the pairwise tail dependence coefficints.
            
    # Example:
    #   LPP = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   setTailRisk <- .nigDependencyFit(LPP)
    
    # Arguments:
    
    # FUNCTION:
    
    # Tail Risk ?
    spec@model$tailRisk = value  
    
    # Return Value:
    spec
}


################################################################################
 

"setSolver<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets the solver value for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
      
    # Set Solver:
    spec@optim$solver = value
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setObjective<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets the solver objective function name for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
      
    # Set Solver:
    spec@optim$objective = value
    
    # Return Value:
    spec
}


# ------------------------------------------------------------------------------


"setTrace<-" <- 
    function(spec, value)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets the trace value for a portfolio structure
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Trace:
    spec@optim$trace = value
    
    # Return Value:
    spec
}


################################################################################

