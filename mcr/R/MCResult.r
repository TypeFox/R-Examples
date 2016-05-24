###############################################################################
##
## MCResult.R
##
## Definition of class MCResult
## Base class of mcreg result objects.
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
    Class="MCResult",
    representation=representation(
        ## measurement data in wide format, one pair of observations per sample:
        ## format: sid, x, y (samples id, method1 measurement, method2 measurement)
        data = "data.frame",
        
        ## regression parameters
        ## Rows: Intercept, Slope
        ## Cols: EST, SE, LCI, UCI
        para = "matrix",

        ## Method annotation
        mnames = "character",

        ## Used regression method
        ## "LinReg","QWLinreg","Deming","PaBa"
        regmeth = "character",

        ## Used methid for CI-calculation
        ## "analytical", "jackknife","bootstrap","nested bootstrap"
        cimeth = "character",
  	
        ## Error Ratio
        error.ratio = "numeric",

        ## Confidence level
        alpha = "numeric",
        
        ## Weights
        weight = "numeric" 
    )
)


###############################################################################
## Method registration
###############################################################################

setMethod(f="initialize",signature="MCResult",definition=MCResult.initialize)

setGeneric("getCoefficients",function(.Object,...){standardGeneric("getCoefficients")})
setMethod("getCoefficients",signature=c(.Object="MCResult"),definition=MCResult.getCoefficients)

setGeneric("getData",function(.Object,...){standardGeneric("getData")})
setMethod("getData",signature=c(.Object="MCResult"),definition=MCResult.getData)

setGeneric("getErrorRatio",function(.Object,...){standardGeneric("getErrorRatio")})
setMethod("getErrorRatio",signature=c(.Object="MCResult"),definition=MCResult.getErrorRatio)

setGeneric("getWeights",function(.Object,...){standardGeneric("getWeights")})
setMethod("getWeights",signature=c(.Object="MCResult"),definition=MCResult.getWeights)

setGeneric("getResiduals",function(.Object,...){standardGeneric("getResiduals")})
setMethod("getResiduals",signature=c(.Object="MCResult"),definition=MCResult.getResiduals)

setGeneric("getFitted",function(.Object,...){standardGeneric("getFitted")})
setMethod("getFitted",signature=c(.Object="MCResult"),definition=MCResult.getFitted)

setGeneric("getRegmethod",function(.Object,...){standardGeneric("getRegmethod")})
setMethod("getRegmethod",signature=c(.Object="MCResult"),definition=MCResult.getRegmethod)

setGeneric("calcCUSUM",function(.Object,...){standardGeneric("calcCUSUM")})
setMethod("calcCUSUM",signature=c(.Object="MCResult"),definition=MCResult.calcCUSUM)

setGeneric("plotDifference",function(.Object,...){standardGeneric("plotDifference")})
setMethod(f="plotDifference",signature=c("MCResult"),definition=MCResult.plotDifference)

setGeneric("calcResponse",function(.Object,...){standardGeneric("calcResponse")})
setMethod(f="calcResponse",signature=c("MCResult"),definition=MCResult.calcResponse)

setGeneric("calcBias",function(.Object,...){standardGeneric("calcBias")})
setMethod("calcBias",signature=c(.Object="MCResult"),definition=MCResult.calcBias)

setMethod("plot",signature=c(x="MCResult"),definition=MCResult.plot)

setGeneric("plotBias",function(x,...){standardGeneric("plotBias")})
setMethod("plotBias",signature=c(x="MCResult"),definition=MCResult.plotBias)

setGeneric("printSummary",function(.Object,...){standardGeneric("printSummary")})
setMethod("printSummary",signature=c(.Object="MCResult"),definition=MCResult.printSummary)

setGeneric("plotResiduals",function(.Object,...){standardGeneric("plotResiduals")})
setMethod("plotResiduals",signature=c(.Object="MCResult"),definition=MCResult.plotResiduals)


