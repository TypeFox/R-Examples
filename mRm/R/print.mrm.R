#    Function: Print routine for object of class mrm.
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

print.mrm <- function(x, ...){
  cat("\n")
  cat("Results - mRm: \n")
  cat("\n")
  cat("Conditional log-likelihood:", x$logLik, "\n")
  cat("BIC:", x$BIC, "\n")
  cat("Number of iterations:", x$number.of.iterations, "\n")
  cat("Number of parameters:", x$number.of.parameters, "\n")
  cat("Convergence to Boundary:", x$conv.to.bound == 1, "\n")
  cat("\n")
  cat("Estimated Item Easiness Parameters:")  
  cat("\n")
  cat("\n")
  para <- x$beta 
  rownames(para) <- paste("Item Nr.", 1:dim(para)[1], sep = " ")
  colnames(para) <- paste("Class Nr.", 1:dim(para)[2], sep = " ")
  print(para)
  cat("\n")
  cat("\n")
  cat("Estimated Class Sizes:")  
  cat("\n")
  cat("\n")
  para <- x$class.size 
  colnames(para) <- paste("Class Nr.", 1:dim(para)[2], sep = " ")
  rownames(para) <- ""
  print(para)
  cat("\n")
  cat("\n")
  cat("Estimated Latent Score Probabilities:")  
  cat("\n")
    cat("\n")
  para <- x$pi.r.c
  rownames(para) <- 1:dim(para)[1]
  colnames(para) <- paste("Class Nr.", 1:dim(para)[2], sep = " ")
  print(para)
 }
