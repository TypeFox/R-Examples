# file NominalLogisticBiplot/R/hermquad.R
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

#Function that calculates the gauss-hermite quadrature in one dimension
#----------------------Parameters--------------
  #N: number of desired points of the quadrature
hermquad <- function(N) {
	EPS = 3e-14
	PIM4 = 0.751125544464943
	MAXIT = 10
	X = matrix(0, N, 1)
	W = matrix(0, N, 1)

	for (i in 1:((N + 1)/2)) {
		if (i == 1)
			z = sqrt(2 * N + 1) - 1.85575 * ((2 * N + 1)^(-0.16667))
		else if (i == 2)
			z = z - (1.14 * N^0.426/z)
		else if (i == 3)
			z = 1.86 * z - 0.86 * X[1]
		else if (i == 4)
			z = 1.91 * z - 0.91 * X[2]
		else z = 2 * z - X[(i - 2)]

		for (iter in 1:(MAXIT + 1)) {
			p1 = PIM4
			p2 = 0
			for (j in 1:N) {
				p3 = p2
				p2 = p1
				p1 = z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
			}
			# the derivative
			pp = sqrt(2 * N) * p2
			# newton step
			z1 = z
			z = z1 - p1/pp
			if (abs(z - z1) <= EPS)
				break
		}
		if (iter == MAXIT + 1)
			print("Too many iterations in hermquad")

		X[i] = z
		X[N + 1 - i] = -z
		W[i] = 2/(pp * pp)
		W[N + 1 - i] = W[i]
	}
	QUAD = list()
	QUAD$X = X
	QUAD$W = W
	class(QUAD)<-"GaussQuadrature"
	return(QUAD)
}
