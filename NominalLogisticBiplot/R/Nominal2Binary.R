# file NominalLogisticBiplot/R/Nominal2Binary.R
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

#This function transfomates a nominal variable into a matrix with so many
#  colums as categories of that variable and each row of the matrix has
#  only a single value one for that level of the category it presents.
#----------------------Parameters--------------
  #y: A nominal variable measured for a set of individuals
Nominal2Binary <- function(y){
	ncat=max(y)
	n=length(y)
	Z=matrix(0,n,ncat)
	for (i in 1:n)
	 Z[i,y[i]]=1

	return(Z)
}
