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

LKrig.cov <- function(x1, x2 = NULL, LKinfo, C = NA, marginal = FALSE) {
	PHI1 <- LKrig.basis(x1, LKinfo)
	# sparse precision matrix for the basis coeficients	
	Q <- LKrig.precision(LKinfo)
	Qc <- chol(Q, memory = LKinfo$choleskyMemory)
	# note: construction of lattice basis depends on alpha and a.wght  and normalizes
	# the basis 
if (!marginal) {
		if (is.null(x2)) {
			PHI2 <- PHI1
		} else {
			PHI2 <- LKrig.basis(x2, LKinfo)
		}
		if (is.na(C[1])) {
			A <- forwardsolve(Qc, transpose = TRUE, t(PHI2), upper.tri = TRUE)
			A <- backsolve(Qc, A)
			return(PHI1 %*% A)
		} else {
			A <- forwardsolve(Qc, transpose = TRUE, t(PHI2) %*% C, upper.tri = TRUE)
			A <- backsolve(Qc, A)
			return(PHI1 %*% A)
		}
	} else {
		if (!is.null(x2)) {
			stop("x2 should not be passed to find marginal variance")
		}
#  NOTE: if LKinfo$normalize = TRUE then basis functions
#   will be already normalized so that the marginal varinace is one
#   without the additional factor of rho.
		PHI <- LKrig.basis(x1, LKinfo)
		marginal.variance <- LKrig.quadraticform(LKrig.precision(LKinfo), 
			PHI, choleskyMemory = LKinfo$choleskyMemory)
		if (!is.null(LKinfo$rho.object)) {
			# add in additional scaling if part of covariance model
			marginal.variance <- marginal.variance * predict(LKinfo$rho.object, 
				x1)
		}
		return(marginal.variance)
	}
	# should not get here.
	}

