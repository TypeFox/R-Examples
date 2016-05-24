###############################################################################
##
## MCResultResampling.r
##
## Definition of class MCResultResampling
## Class of mcreg result objects that contain Bootstrap based results.
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
  Class = "MCResultResampling",
  representation = representation(
      # Regression coefficients for global data
      # Intersept
      # Slope
      glob.coef = "numeric",
      glob.sigma = "numeric",
      
      ## global mean (weighted mean)
      xmean = "numeric",

      # Number of samples for Bootstrap,
      # for method = "jackknife" it as always n (length of dataset)
      nsamples = "numeric",

      # Number of nested samples for nested Bootstrap
      nnested = "numeric",

      # Intercept for each bootstrap sample
      # length of B0 must be equal nsamples
      B0 = "numeric",

      # Slope for each bootstrap sample,
      # length of B1 must be equal nsamples
      B1 = "numeric",

      # Standard deviation of sample coefficients
      # (for nested bootstrap it is a vector)
      sigmaB0 = "numeric",
      sigmaB1 = "numeric",
      
      # Mean X for bootstrap sample      
      MX = "numeric",

      # bootstrap method for CI calculation
      bootcimeth="character",
      
      # RNG settings
      rng.seed="numeric",
      rng.kind="character"
  ),
   
   contains = "MCResult"
)


###############################################################################
## Method registration
###############################################################################

setMethod(f="initialize",signature="MCResultResampling",definition=MCResultResampling.initialize)

setGeneric("plotBootstrapCoefficients",function(.Object,...){standardGeneric("plotBootstrapCoefficients")})
setMethod("plotBootstrapCoefficients",signature=c(.Object="MCResultResampling"),definition=MCResultResampling.plotBootstrapCoefficients)

setGeneric("plotBootstrapT",function(.Object,...){standardGeneric("plotBootstrapT")})
setMethod("plotBootstrapT",signature=c(.Object="MCResultResampling"),definition=MCResultResampling.plotBootstrapT)

setGeneric("bootstrapSummary",function(.Object,...){standardGeneric("bootstrapSummary")})
setMethod("bootstrapSummary",signature=c(.Object="MCResultResampling"),definition=MCResultResampling.bootstrapSummary)

setMethod("calcResponse",signature=c(.Object="MCResultResampling"),definition=MCResultResampling.calcResponse)

setMethod(f="printSummary",signature=c(.Object="MCResultResampling"),definition=MCResultResampling.printSummary)

