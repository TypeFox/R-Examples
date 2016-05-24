
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
# FUNCTION:                DESCRIPTION:
#  getStatistics            Extracts statistics slot
#   getMean                  Extracs mean from statistics
#   getCov                   Extracs covariance Sigma from statistics
#   getMu                    Extracs mu from statistics
#   getSigma                 Extracs Sigma from statistics
#   getEstimator             Extracts estimator from 
# FUNCTION:                DESCRIPTION:
#  getTailRisk               Extracts tailRisk slot
################################################################################


# fPFOLIODATA:

# data = list(
#   series
#   nAssets
#   names)

# statistics = list(
#   mean,
#   Cov,
#   mu,
#   Sigma,
#   estimator) 

# tailRisk = list()  


# ------------------------------------------------------------------------------t


getData.fPFOLIODATA <- function(object) object@data
# Extracts the @data slot from a fPFOLIODATA object
getSeries.fPFOLIODATA <- function(object) object@data$series   
getNAssets.fPFOLIODATA <- function(object) object@data$nAssets
getUnits.fPFOLIODATA <- function(x) x@data$names


# ------------------------------------------------------------------------------


getStatistics.fPFOLIODATA <- function(object) object@statistics
# Extracts the @statistics slot from a fPFOLIODATA object
getMean.fPFOLIODATA <- function(object) object@statistics$mean
getCov.fPFOLIODATA <- function(object) object@statistics$Cov
getEstimator.fPFOLIODATA <- function(object) object@statistics$estimator
getMu.fPFOLIODATA <- function(object) object@statistics$mu
getSigma.fPFOLIODATA <- function(object) object@statistics$Sigma


# ------------------------------------------------------------------------------


getTailRisk.fPFOLIODATA <- function(object) object@tailRisk
# Extracts the @tailRisk slot from a fPFOLIODATA object 

################################################################################

