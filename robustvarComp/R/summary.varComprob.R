#############################################################
# 
#	summary functions
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: July, 01, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

summary.varComprob <- function(object, print.outliers=FALSE, ...) {
  beta <- object$beta
  se.beta <- sqrt(diag(object$vcov.beta))
  z.beta <- beta/se.beta
  p.value <- 2*pnorm(-abs(z.beta))
  object$zTable <- cbind(beta, se.beta, z.beta, p.value)
  colnames(object$zTable) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object <- c(object, summarizeRobMahalanobis(object, ...))
  object$print.outliers <- print.outliers
  class(object) <- c("summary.varComprob", class(object))
  return(object)
}

summary.varComprob.compositeS <- summary.varComprob.compositeTau <- summary.varComprob.compositeMM <- summary.varComprob.S <- summary.varComprob.Tau <- summary.varComprob.MM <- summary.varComprob

print.summary.varComprob <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), print.outliers=x$print.outliers, ...) {
  cat("Method: ", x$control$method, "\n\n")
  cat("Fixed effects: \n")
  printCoefmat(x$zTable, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  cat("\n Random effect variances: \n")
  print.default(x$eta)
  cat("\n Residual variance: \n")
  print.default(x$eta0)  
  cat("\n Value of the objective function at convergence: \n")
  print.default(x$min)

  if (is.na(pmatch(x="composite", table=x$control$method))) {  
    cat("\n")
    cat("\n")
    summarizeRobWeights(x$weights)
  }
  
  if (print.outliers) {
    cat("\n")
    cat("\n")
  
    cat("Cell outliers:\n")
    if (any(x$cellsoutliers$pos.cells)) {
      print.default(x$cellsoutliers$table.cells, digits=digits, ...)     
      cat(" on", x$cellsoutliers$num.cells, "cells\n")
    } else
      cat("No cell outliers were found\n")
    cat("\n")
    cat("Couple outliers:\n")
    if (any(x$couplesoutliers$pos.couples)) {
      print.default(x$couplesoutliers$table.couples, digits=digits, ...)     
      cat(" on", x$couplesoutliers$num.couples, "couples\n")
    } else
      cat("No couple outliers were found\n")
    cat("\n")
    cat("Row outliers\n")
    if (any(x$rowsoutliers$pos.rows)) {
      print.default(x$rowsoutliers$table.rows, digits=digits, ...)     
      cat(" on", x$rowsoutliers$num.rows, "observations\n")
    } else
      cat("No row outliers were found\n")
  }
  invisible(x)
}

print.summary.varComprob.compositeS <- print.summary.varComprob.compositeTau <- print.summary.varComprob.compositeMM <- print.summary.varComprob.S <- print.summary.varComprob.Tau <- print.summary.varComprob.MM <- print.summary.varComprob
