
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA 02111-1307 USA


################################################################################
# FUNCTION:                DESCRIPTION:                                         
#  getData                  Extracts data slot                                  
#   getSeries                Extracts assets series data                        
#   getNAssets               Extracts number of assets from data                
#   getUnits                 Extracts assets names from data                    
#  getStatistics            Extracts statistics slot                            
#   getMean                  Extracs mean from statistics                       
#   getCov                   Extracs covariance Sigma from statistics           
#   getMu                    Extracs mu from statistics                         
#   getSigma                 Extracs Sigma from statistics                      
#   getEstimator             Extracts estimator from                            
#  getTailRisk              Extracts tailRisk slot                              
# FUNCTION:                DESCRIPTION:                                         
#  getSpec                  Extracs specification Slot                          
#   getType                  Extracts type of portfolio                         
#   getOptimize              Extracts what to optimize of portfolio             
#   getEstimator             Extracts mean-covariance estimator                 
#   getParams                Extracts optional parameter list                   
#    getAlpha                 Extracts target VaR-alpha specification           
#    getA                     Extracts quadratic LPM exponent specification     
#  getPortfolio             Extract portfolio slot                              
#   getWeights               Extracts weights from a portfolio object           
#   getTargetReturn          Extracts target return from specification          
#   getTargetRisk            Extracts target riks from specification            
#   getRiskFreeRate          Extracts risk free rate from specification         
#   getNFrontierPoints       Extracts number of frontier points                 
#   getStatus                Extracts portfolio status information              
#  getOptim                 Extract optim slot                                  
#   getSolver                Extracts solver from specification                 
#   getObjective             Extracts objective
#   getOptions               Extracts optimization options
#   getControl               Extracts solver control options
#   getTrace                 Extracts solver's trace flag
# FUNCTION:                 DESCRIPTION:
#  getConstraints            Extracts weight constraints
# FUNCTION:                 DESCRIPTION:
#  getCovRiskBudgets         Extracts covariance risk budgets
#  getTailRiskBudgets        Extracts tail risk budgets
################################################################################


# Extract from data slot of an object of class fPORTFOLIO:

        
getData.fPORTFOLIO <- 
    function(object) object@data

getSeries.fPORTFOLIO <- 
    function(object) object@data@data$series
getNAssets.fPORTFOLIO <- 
    function(object) object@data@data$nAssets
getUnits.fPORTFOLIO <- 
  function(x) x@data@data$names
    
    
getStatistics.fPORTFOLIO <- 
    function(object) object@data@statistics 
getMean.fPORTFOLIO <- 
    function(object) object@data@statistics$mean
getCov.fPORTFOLIO <- 
    function(object) object@data@statistics$Cov
getEstimator.fPORTFOLIO <- 
    function(object) object@data@statistics$estimator
getMu.fPORTFOLIO <- 
    function(object) object@data@statistics$mu
getSigma.fPORTFOLIO <- 
    function(object) object@data@statistics$Sigma
 
    
# ------------------------------------------------------------------------------


# Extract from spec slot of an object of class fPORTFOLIO:


getSpec.fPORTFOLIO <- 
    function(object) object@spec 
getModel.fPORTFOLIO <- 
    function(object) object@spec@model
getType.fPORTFOLIO <- 
    function(object) object@spec@model$type
getOptimize.fPORTFOLIO <- 
    function(object) object@spec@model$optimize
getEstimator.fPORTFOLIO <- 
    function(object) object@spec@model$estimator
getTailRisk.fPORTFOLIO <- 
    function(object) object@spec@model$tailRisk
getParams.fPORTFOLIO <- 
    function(object) object@spec@model$params
getAlpha.fPORTFOLIO <- 
    function(object) object@spec@model$params$alpha
getA.fPORTFOLIO <- 
    function(object) object@spec@model$params$a


# DW object@spec renamed to object@portfolio
getPortfolio.fPORTFOLIO <-     
    function(object) object@portfolio@portfolio
getWeights.fPORTFOLIO <-      
    function(object) object@portfolio@portfolio$weights  
getTargetReturn.fPORTFOLIO <- 
    function(object) object@portfolio@portfolio$targetReturn
getTargetRisk.fPORTFOLIO <-   
    function(object) object@portfolio@portfolio$targetRisk
getRiskFreeRate.fPORTFOLIO <- 
    function(object) object@spec@portfolio$riskFreeRate
getNFrontierPoints.fPORTFOLIO <- 
    function(object) object@spec@portfolio$nFrontierPoints
getStatus.fPORTFOLIO <-  
    function(object) object@spec@portfolio$status

    
getOptim.fPORTFOLIO <-      
    function(object) object@spec@optim
getSolver.fPORTFOLIO <-    
    function(object) object@spec@optim$solver
getObjective.fPORTFOLIO <- 
    function(object) object@spec@optim$objective 
getOptions.fPORTFOLIO <-   
    function(object) object@spec@optim$options 
getControl.fPORTFOLIO <-   
    function(object) object@spec@optim$control 
getTrace.fPORTFOLIO <-     
    function(object) object@spec@optim$trace

getCovRiskBudgets.fPORTFOLIO <- 
    function(object) object@portfolio@portfolio$covRiskBudgets


# ------------------------------------------------------------------------------


# Extract from constraints slot of an object of class fPORTFOLIO:


getConstraints.fPORTFOLIO <- 
    function(object) object@constraints@stringConstraints
    
    
getConstraintsTypes <- 
    function(object) 
{
    Constraints = getConstraints(object)
    Types = NULL
    if(!is.na(pmatch("LongOnly", Constraints))) Types = c(Types, "LongOnly") 
    if(!is.na(pmatch("Short", Constraints))) Types = c(Types, "Short") 
    if(!is.na(pmatch("minW", Constraints))) Types = c(Types, "minW") 
    if(!is.na(pmatch("maxW", Constraints))) Types = c(Types, "maxW") 
    if(!is.na(pmatch("minsumW", Constraints))) Types = c(Types, "minsumW") 
    if(!is.na(pmatch("maxsumW", Constraints))) Types = c(Types, "maxsumW") 
    if(!is.na(pmatch("minB", Constraints))) Types = c(Types, "minB") 
    if(!is.na(pmatch("maxB", Constraints))) Types = c(Types, "maxB") 
    Types
}


################################################################################


.getCovRiskBudgets.fPORTFOLIO <-  
function (object) 
{   # A function implemented by Rmetrics

    # Description:
    #   Extracts risk budgets from a portfolio object
    
    # FUNCTION:
    
    # Covariance Risk Budgets:
    weights = object@portfolio$weights
    ans = NA
    Sigma = object@data$data@statistics$Sigma
    if (is.null(dim(weights))) {
        # Single Portfolio ...
        ans1 = as.vector(weights %*% Sigma %*% weights)
        ans2 = as.vector(weights * Sigma %*% weights)
        ans = round(ans2/ans1, digits = 4)
        names(ans) = names(weights)
    } else {
        # Frontier ...
        Names = colnames(weights)
        ans = NULL
        for (i in 1:(dim(weights)[1])) {
            ans1 = as.vector(weights[i, ] %*% Sigma %*% weights[i, ])
            ans2 = as.vector(weights[i, ] * Sigma %*% weights[i, ])
            ans = rbind(ans, ans2/ans1)
        }
        colnames(ans) = Names
    }
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


getTailRiskBudgets.fPORTFOLIO <-  
function (object) 
{   # A function implemented by Rmetrics

    # Description:
    #   Extracts tail risk budgets from a portfolio object
    
    # Arguments:
    #   object - an object of S4 class fPORTFOLIO as returned by the
    #       functions *Portfolio().
    
    # FUNCTION:
    
    # Check if available:
    Lambda = object@spec@model$tailRisk$lower
    if (is.null(Lambda)) return(NA)
    
    # Tail Risk Budgets:
    weights = getWeights(object)
    ans = NA
    if (is.null(dim(weights))) {
        ans1 = as.vector(weights %*% Lambda %*% weights)
        ans2 = as.vector(weights * Lambda %*% weights)
        ans1 = 1
        ans = round(ans2/ans1, digits = 4)
        names(ans) = names(weights)
    }
    else {
        Names = colnames(weights)
        ans = NULL
        for (i in 1:(dim(weights)[1])) {
            ans1 = as.vector(weights[i, ] %*% Lambda %*% weights[i, ])
            ans2 = as.vector(weights[i, ] * Lambda %*% weights[i, ])
            ans1 = 1
            ans = rbind(ans, ans2/ans1)
        }
        colnames(ans) = Names
    }
    
    # Return Value:
    ans
}


################################################################################

