# file NominalLogisticBiplot/R/NominalDistances.R
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

#Function that calculates a distance matrix starting from a data set with some nominal variables.
#----------------------Parameters--------------
  #x: data matrix with the information of the nominal variables.
  #similarities: boolean parameter to choose if we want the similarities or the distances.
NominalDistances <- function(x, similarities = FALSE) {
	n <- nrow(x)
	p <- ncol(x)
	sim = matrix(1, n, n)
	distance = matrix(0, n, n)
	for (i in 1:n) {
		for (j in i:n) {
      #print(c(i,j))
			sim[i, j] = sum(as.numeric(x[i, ] == x[j, ]))/p
			sim[j, i] = sim[i, j]
		}
	}

	distance=sqrt(1-sim)
	if (similarities)
		return(sim)
	else return(distance)
}


