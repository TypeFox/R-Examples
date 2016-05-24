
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
# FUNCTION:                     DESCRIPTION:
#  getA                          Defines Use Method for A                
#  getAlpha                      Defines Use Method for Alpha            
#  getConstraints                Defines Use Method for Constraints      
#  getControl                    Defines Use Method for Control          
#  getCov                        Defines Use Method for Cov              
#  getCovRiskBudgets             Defines Use Method for CovRiskBudgets   
#  getData                       Defines Use Method for Data             
#  getEstimator                  Defines Use Method for Estimator        
#  getMean                       Defines Use Method for Mean             
#  getMu                         Defines Use Method for Mu               
#  getNAssets                    Defines Use Method for NAssets          
# .getNames                      Defines Use Method for Names            
#  getNFrontierPoints            Defines Use Method for NFrontierPoints  
#  getMessages                   Defines Use Method for Messages         
#  getObjective                  Defines Use Method for Objective        
#  getOptim                      Defines Use Method for Optim            
#  getOptimize                   Defines Use Method for Optimize         
#  getOptions                    Defines Use Method for Options          
#  getPortfolio                  Defines Use Method for Portfolio        
#  getParams                     Defines Use Method for Params           
#  getRiskFreeRates              Defines Use Method for RiskFreeRates    
#  getSeries                     Defines Use Method for Series           
#  getSigma                      Defines Use Method for Sigma            
#  getSolver                     Defines Use Method for Solver           
#  getSpec                       Defines Use Method for Spec             
#  getStatistics                 Defines Use Method for Statistics       
#  getStatus                     Defines Use Method for Status           
#  getTailRisk                   Defines Use Method for TailRisk         
#  getTailRiskBudgets            Defines Use Method for TailRiskBudgets  
#  getTargetReturn               Defines Use Method for TargetReturn     
#  getTargetRisk                 Defines Use Method for TargetRisk       
#  getTrace                      Defines Use Method for Trace            
#  getType                       Defines Use Method for Type  
#  getUnits                      Defines Use Method for Units [Asset Names]           
#  getWeights                    Defines Use Method for Weights          
################################################################################


getA <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getA")
}


# ------------------------------------------------------------------------------


getAlpha <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getAlpha")
}


# ------------------------------------------------------------------------------


getConstraints <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getConstraints")
}


# ------------------------------------------------------------------------------


getControl <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getControl")
}


# ------------------------------------------------------------------------------


getCov <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getCov")
}


# ------------------------------------------------------------------------------


getData <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getData")
}


# ------------------------------------------------------------------------------


getCovRiskBudgets <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getCovRiskBudgets")
}


# ------------------------------------------------------------------------------


getEstimator <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getEstimator")
}


# ------------------------------------------------------------------------------


getMean <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Return Value:
    UseMethod("getMean")
}


# ------------------------------------------------------------------------------

getMu <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getMu")
}


# ------------------------------------------------------------------------------


getNAssets <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Return Value:
    UseMethod("getNAssets")
}


# ------------------------------------------------------------------------------


.getNames <-
    function(object)
{
    # A function implemented by Diethelm Wuertz
    
    # REPLACED BY getUnits.

    # Return Value:
    UseMethod("getNames")
}


# ------------------------------------------------------------------------------


getNFrontierPoints <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getNFrontierPoints")
}


# ------------------------------------------------------------------------------


getMessages <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getMessages")
}


# ------------------------------------------------------------------------------


getObjective <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    UseMethod("getObjective")
}


# ------------------------------------------------------------------------------


getOptim <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getOptim")
}


# ------------------------------------------------------------------------------


getOptimize <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    UseMethod("getOptimize")
}


# ------------------------------------------------------------------------------


getOptions <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getOptions")
}


# ------------------------------------------------------------------------------


getPortfolio <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Return Value:
    UseMethod("getPortfolio")
}


# ------------------------------------------------------------------------------


getParams <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    UseMethod("getParams")
}


# ------------------------------------------------------------------------------


getRiskFreeRate <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getRiskFreeRate")
}


# ------------------------------------------------------------------------------
# DW: Take care of getSeries in package timeSeries


getSeries <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getSeries")
}


# ------------------------------------------------------------------------------


getSigma <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getSigma")
}


# ------------------------------------------------------------------------------


getSolver <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getSolver")
}


# ------------------------------------------------------------------------------


getSpec <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getSpec")
}


# ------------------------------------------------------------------------------


getStatistics <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getStatistics")
}


# ------------------------------------------------------------------------------


getStatus <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getStatus")
}


# ------------------------------------------------------------------------------


getTailRisk <-
    function(object)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    UseMethod("getTailRisk")
}


# ------------------------------------------------------------------------------


getTailRiskBudgets <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getTailRiskBudgets")
}


# ------------------------------------------------------------------------------


getTargetReturn <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getTargetReturn")
}


# ------------------------------------------------------------------------------


getTargetRisk <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getTargetRisk")
}


# ------------------------------------------------------------------------------


getTrace <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getTrace")
}


# ------------------------------------------------------------------------------


getType <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getType")
}


# ------------------------------------------------------------------------------
# DW: already defined in package timeSeries

# getUnits <-
#     function(object)
# {
#     # A function implemented by Diethelm Wuertz
# 
#     # FUNCTION: 
#     
#     # Return Value:
#     UseMethod("getUnits")
# }


# ------------------------------------------------------------------------------
# TS: already defined in package fBasics

# getModel <-
#   function(object)
#   {
#     # A function implemented by Tobias Setz
#     
#     # FUNCTION: 
#     
#     # Return Value:
#     UseMethod("getModel")
#   }


# ------------------------------------------------------------------------------


getWeights <-
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION: 
    
    # Return Value:
    UseMethod("getWeights")
}


################################################################################

