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

LKrig.sim <- function(x1, LKinfo, M = 1, just.coefficients = FALSE) {
	Q <- LKrig.precision(LKinfo)
	#
	#  Q is precision matrix of the coefficients -- not of the field
#  last step does the multiplication to go from coefficients to evaluating
#  values at the field
#  Q = t(H)%*%H = inv((Sigma)
#  So   Sigma= Hi%*% t(Hi)
#  find u= t(Hi) %*% N(0,1)   then cov(u) = t(Hi)%*%Hi
#  Hi is upper triangular
#
# snippet of code to test the algebra ...
#   x<-seq( 0,1,,20); Sigma<- exp(-rdist( x,x)/2.5); Q<- solve( Sigma)
#   Mc <- chol(Q); H<- Mc ; Hi<- solve(H);
#   test.for.zero( Q, t(H)%*%H); test.for.zero(Sigma, Hi%*%t(Hi))
#   E<- rnorm(20);  u1<- Hi%*% E ;   u2<-backsolve(Mc,E)
#   test.for.zero(u1,u2)
#
   Qc <- chol(Q, memory = LKinfo$choleskyMemory)
	m <- LKinfo$latticeInfo$m
	E <- matrix(rnorm(M * m), nrow = m, ncol = M)
	A <- backsolve(Qc, E)
	if (just.coefficients) {
		return(A)
	} else {
		PHI1 <- LKrig.basis(x1, LKinfo)
		return(PHI1 %*% A)
	}
}

