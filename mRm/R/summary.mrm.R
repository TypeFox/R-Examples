#    Function: Summary function for object of class 'mrm'.
#    Copyright (C) 2011  David Preinerstorfer
#    david.preinerstorfer@univie.ac.at
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

`summary.mrm` <- function(object, ...){
  cat("\n")
  cat("Summary of Results - mRm: \n")
  cat("\n")
  cat("Conditional log-likelihood:", object$logLik, "\n")
  cat("BIC:", object$BIC, "\n")
  cat("Number of iterations:", object$number.of.iterations, "\n")
  cat("Number of parameters:", object$number.of.parameters, "\n")
  cat("Number of items:", dim(object$beta)[1], "\n")
  cat("Number of classes:", dim(object$beta)[2], "\n")
  cat("Convergence to Boundary:", object$conv.to.bound == 1, "\n")
  cat("\n")
 }
