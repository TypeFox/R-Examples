# file NominalLogisticBiplot/R/PCoA.R
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

#This function calculates principal coordinates analysis for a distance matrix.
#----------------------Parameters--------------
  #dis: distance matrix
  #r : number of dimensions retained
PCoA <- function(dis, r = 2) {
	n <- dim(dis)[1]
	b <- -0.5 * (diag(n) - matrix(1, n, n)/n) %*% dis^2 %*% (diag(n) - matrix(1, n, n)/n)
	solut <- svd(b)
	Inertia = (solut$d/sum(solut$d)) * 100
	g <- solut$u[1:n, 1:r] %*% diag(sqrt(solut$d[1:r]))
	st <- rowSums(t(solut$u %*% diag(sqrt(solut$d)))^2)
	qlr <- diag(1/st) %*% (solut$u %*% diag(sqrt(solut$d)))^2
	result = list(EigenValues = solut$d, Inertia = Inertia, RowCoordinates = g, RowQualities = qlr[1:n, 1:r])
	return(result)
}
