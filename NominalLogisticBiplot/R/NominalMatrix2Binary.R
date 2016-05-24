# file NominalLogisticBiplot/R/NominalMatrix2Binary.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
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

#This function transfomates a matrix with nominal variables into a matrix
#  with so many colums as the sum for all the level categories of the whole
#  variable data set and each row of the matrix has so many one values
#  as number of variables, situated in the level categories it presents.
#----------------------Parameters--------------
  #Y: A matrix with nominal variables measured for a set of individuals.
NominalMatrix2Binary <- function(Y){
	n=dim(Y)[1]
	p=dim(Y)[2]
	G= Nominal2Binary(Y[,1])
	for (j in 2:p)
	 G=cbind(G, Nominal2Binary(Y[,j]))

	return(G)
}

