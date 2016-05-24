# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

print.LatticeKrig <- function(x, digits = 4, ...) {
	LKinfo <- x$LKinfo
	if (is.matrix(x$residuals)) {
		n <- nrow(x$residuals)
		NData <- ncol(x$residuals)
	} else {
		n <- length(x$residuals)
		NData <- 1
	}
	c1 <- "Number of Observations:"
	c2 <- n
	if (NData > 1) {
		c1 <- c(c1, "Number of data sets fit:")
		c2 <- c(c2, NData)
	}
	c1 <- c(c1, "Number of parameters in the fixed component")
	c2 <- c(c2, x$nt)
	if (x$nZ > 0) {
		c1 <- c(c1, "Number of covariates")
		c2 <- c(c2, x$nZ)
	}
	if (!is.na(x$eff.df)) {
		c1 <- c(c1, " Effective degrees of freedom (EDF)")
		c2 <- c(c2, signif(x$eff.df, digits))
		c1 <- c(c1, "   Standard Error of EDF estimate: ")
		c2 <- c(c2, signif(x$trA.SE, digits))
	}
	c1 <- c(c1, "MLE lambda = sigma^2/rho ")
	c2 <- c(c2, signif(x$lambda, digits))
	if (NData == 1) {
		c1 <- c(c1, "MLE sigma ")
		c2 <- c(c2, signif(x$shat.MLE, digits))
		c1 <- c(c1, "MLE rho")
		c2 <- c(c2, signif(x$rho.MLE, digits))
	}

	sum <- cbind(c1, c2)
	dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
	cat("Call:\n")
	dput(x$call)
	print(sum, quote = FALSE)
	cat(" ", fill = TRUE)
	cat("sigma and rho found by maximum likelihood for ", fill = TRUE)
	cat("fixed values of alpha and a.wght", fill = TRUE)
	cat(" ", fill = TRUE)
	#    
	if (is.null(x$LKinfo$fixedFunction)) {
		cat("No fixed part of model", fill = TRUE)
	} else {

		if (x$LKinfo$fixedFunction == "LKrigDefaultFixedFunction") {
			cat("Fixed part of model is a polynomial of degree", x$LKinfo$fixedFunctionArgs$m - 
				1, "(m-1)", fill = TRUE)
		} else {
			cat("Fixed part of model uses the function:", x$LKinfo$fixedFunction, 
				fill = TRUE)
			cat("with the argument list:", fill = TRUE)
			print(x$LKinfo$fixedFunctionArgs)
		}
	}
	cat("Basis function : ", LKinfo$basisInfo$BasisType, 
			fill = TRUE)
	cat("Basis function used: ", LKinfo$basisInfo$BasisFunction, 
			fill = TRUE)
	cat("Distance metric: ", LKinfo$distance.type, fill = TRUE)
	cat(" ", fill = TRUE)
	cat("Lattice summary:", fill = TRUE)
	cat(LKinfo$nlevel, " Level(s)", LKinfo$latticeInfo$m, "basis functions", 
		"with overlap of ", LKinfo$basisInfo$overlap, "(lattice units)", 
		fill = TRUE)
	cat(" ", fill = TRUE)
	#
	temp <- cbind(1:LKinfo$nlevel, LKinfo$latticeInfo$mLevel, LKinfo$latticeInfo$delta)
	dimnames(temp) <- list(rep("", LKinfo$nlevel), c("Level", "Lattice points", 
		"Spacing"))
	print(temp)
	cat(" ", fill = TRUE)
	cat("Nonzero entries in Ridge regression matrix", x$nonzero.entries, 
		fill = TRUE)	
	#
	cat(" ", fill = TRUE)
	if (length(LKinfo$alpha[[1]]) == 1) {
		cat("Value(s) for alpha: ", unlist(LKinfo$alpha), fill = TRUE)
	} else {
		cat("alpha values passed as a vector for each level", fill = TRUE)
	}
	#    
	cat(" ", fill = TRUE)
	if (length(LKinfo$a.wght[[1]]) == 1) {
		a.wght <- unlist(LKinfo$a.wght)
		cat("Value(s) for lattice dependence (a.wght): ", a.wght, fill = TRUE)
	} else {
		cat("Value(s) for weighting in GMRF (a.wght): ", unlist(LKinfo$alpha), 
			fill = TRUE)
	}
	#    
	cat(" ", fill = TRUE)
	if (LKinfo$normalize) {
		cat("Basis functions normalized so marginal process variance is stationary", 
			fill = TRUE)
	}
	invisible(x)
}

