###############################################################################
##
## MCResultJackknife.R
##
## Definition of class MCResultJackknife
## Class of mcreg result objects that contain Jackknife based results.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

###############################################################################
## Class definition
###############################################################################

setClass(
    Class = "MCResultJackknife",
    representation = representation(
        # Regression coefficients for global data
        # Intersept
        # Slope
        glob.coef = "numeric",
        
        B0jack = "numeric",
        B1jack = "numeric"
    ),
    contains = "MCResult"
)


###############################################################################
## Method registration
###############################################################################

setMethod(f="initialize",signature="MCResultJackknife",definition=MCResultJackknife.initialize)

setMethod("calcResponse",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.calcResponse)

setGeneric("getJackknifeSlope",function(.Object,...){standardGeneric("getJackknifeSlope")})
setMethod("getJackknifeSlope",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.getJackknifeSlope)

setGeneric("getJackknifeIntercept",function(.Object,...){standardGeneric("getJackknifeIntercept")})
setMethod("getJackknifeIntercept",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.getJackknifeIntercept)

setGeneric("getRJIF",function(.Object,...){standardGeneric("getRJIF")})
setMethod("getRJIF",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.getRJIF)

setGeneric("plotwithRJIF",function(.Object,...){standardGeneric("plotwithRJIF")})
setMethod("plotwithRJIF",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.plotwithRJIF)

setGeneric("getJackknifeStatistics",function(.Object,...){standardGeneric("getJackknifeStatistics")})
setMethod("getJackknifeStatistics",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.getJackknifeStatistics)

setMethod(f="printSummary",signature=c(.Object="MCResultJackknife"),definition=MCResultJackknife.printSummary)

