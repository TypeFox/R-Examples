###############################################################################
##
## MCResultBCa.R
##
## Definition of class MCResultBCa
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
  Class = "MCResultBCa",
  contains = c("MCResultJackknife"),
  representation = representation(
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

      ## ci method, should be BCa
      bootcimeth="character",
      
      # RNG settings
      rng.seed="numeric",
      rng.kind="character"

  )
)


###############################################################################
## Method registration
###############################################################################

setMethod(f="initialize",signature="MCResultBCa",definition=MCResultBCa.initialize)

setMethod("plotBootstrapCoefficients",signature=c(.Object="MCResultBCa"),definition=MCResultBCa.plotBootstrapCoefficients)

setMethod("plotBootstrapT",signature=c(.Object="MCResultBCa"),definition=MCResultBCa.plotBootstrapT)

setMethod("bootstrapSummary",signature=c(.Object="MCResultBCa"),definition=MCResultBCa.bootstrapSummary)

setMethod("calcResponse",signature=c(.Object="MCResultBCa"),definition=MCResultBCa.calcResponse)

setMethod(f="printSummary",signature=c(.Object="MCResultBCa"),definition=MCResultBCa.printSummary)



















