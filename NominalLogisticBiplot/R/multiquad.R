# file NominalLogisticBiplot/R/multiquad.R
# copyright (C) 2012-2013 J.L. Vicente-Villardon and J.C. Hernandez
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#Function that calculates the multidimensional gauss-hermite quadrature.
#----------------------Parameters--------------
  #nnodos: number of nodes for the quadrature for each of the dimensions of the solution
  #dims: number of the dimensions in the rediced space.
multiquad <- function(nnodos, dims) {
	Q = hermquad(nnodos)
	I = patterns_eq(nnodos, dims)
	n = dim(I)[1]
	X2 = matrix(0, n, dims)
	A2 = matrix(1, n, 1)
	for (i in 1:n)
   for (j in 1:dims) {
		X2[i, j] = Q$X[I[i, j]]
		A2[i, 1] = A2[i, 1] * Q$W[I[i, j]]
	 }
	Max = Q$W[1] * Q$W[round((nnodos + 1)/2)]/15
	keep = (A2 > Max)
	n2 = sum(keep)
	X = matrix(0, n2, dims)
	A = matrix(1, n2, 1)
	k = 0
	for (i in 1:n)
   if (keep[i]) {
	   	k = k + 1
		  X[k, ] = X2[i, ]
		  A[k, ] = A2[i, ]
	 }
	QUAD = list()
	QUAD$X = X
	QUAD$A = A
	class(QUAD) <- "MultiGaussQuadrature"
	return(QUAD)
}