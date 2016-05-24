
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
#  getModel                      Extract whole model slot
#   getType                       Extract portfolio type from specification 
#   getOptimize                   Extract what to optimize from specification
#   getEstimator                  Extract type of covariance estimator
#   getTailRisk                   Extract list of tail dependency risk matrixes
#   getParams                     Extract parameters from specification
#    getAlpha                      Extracts target VaR-alpha specification
#    getA                          Extracts quadratic LPM Exponent
# FUNCTION:                     DESCRIPTION:
#  getPortfolio                  Extract whole portfolio slot
#   getWeights                    Extracts weights from a portfolio object
#   getTargetReturn               Extracts target return from specification
#   getTargetRisk                 Extracts target riks from specification
#   getRiskFreeRate               Extracts risk free rate from specification 
#   getNFrontierPoints            Extracts number of frontier points 
#   getStatus                     Extracts portfolio status information
# FUNCTION:                     DESCRIPTION:
#  getOptim                       Extract whole optim slot
#   getSolver                     Extracts solver from specification
#   getObjective                  Extracs name of objective function
#   getOptions                    Extracs options              
#   getControl                    Extracs control list parameters
#   getTrace                      Extracts solver's trace flag
# FUNCTION:                     DESCRIPTION:
#  getMessages                    Extract whole messages slot
################################################################################


# fPFOLIOSPEC:

# model = list(
#   type = "MV",
#   optimize = "minRisk",
#   estimator = "covEstimator",
#   tailRisk = NULL,
#   params = list(alpha = 0.05, a = 1))

# portfolio = list(
#   weights = NULL, 
#   targetReturn = NULL, 
#   targetRisk = NULL, 
#   targetAlpha = NULL,
#   riskFreeRate = 0, 
#   nFrontierPoints = 50,
#   status = 0)

# optim = list(
#   solver = "solveRquadprog",
#   objective = NULL,
#   options = list(meq=2), 
#   control = list(),
#   trace = FALSE)

# messages = list(NULL) 


# ------------------------------------------------------------------------------


getModel.fPFOLIOSPEC <- function(object) object@model
getType.fPFOLIOSPEC <- function(object) object@model$type[1]
getOptimize.fPFOLIOSPEC <- function(object) object@model$optimize
getEstimator.fPFOLIOSPEC <- function(object) object@model$estimator
getTailRisk.fPFOLIOSPEC <- function(object) object@model$tailRisk
getParams.fPFOLIOSPEC <- function(object) object@model$params
getAlpha.fPFOLIOSPEC <- function(object) object@model$params$alpha
getA.fPFOLIOSPEC <- function(object) object@model$params$a 

.getEstimatorFun <- function(object) match.fun(getEstimator(object))


# ------------------------------------------------------------------------------


getPortfolio.fPFOLIOSPEC <- function(object) object@portfolio
getWeights.fPFOLIOSPEC <- function(object) object@portfolio$weights
getTargetReturn.fPFOLIOSPEC <- function(object) object@portfolio$targetReturn
getTargetRisk.fPFOLIOSPEC <- function(object) object@portfolio$targetRisk
getRiskFreeRate.fPFOLIOSPEC <- function(object) object@portfolio$riskFreeRate
getNFrontierPoints.fPFOLIOSPEC <- function(object) object@portfolio$nFrontierPoints
getStatus.fPFOLIOSPEC <-  function(object) object@portfolio$status


# ------------------------------------------------------------------------------


getOptim.fPFOLIOSPEC <- function(object) object@optim
getSolver.fPFOLIOSPEC <- function(object) object@optim$solver 
getObjective.fPFOLIOSPEC <- function(object) object@optim$objective 
getOptions.fPFOLIOSPEC <- function(object) object@optim$options  
getControl.fPFOLIOSPEC <- function(object) object@optim$control
getTrace.fPFOLIOSPEC <- function(object) object@optim$trace


# ------------------------------------------------------------------------------


getMessages.fPFOLIOSPEC <- function(object) object@messages


################################################################################

