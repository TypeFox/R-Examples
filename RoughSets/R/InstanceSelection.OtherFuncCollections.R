#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
# This function is to transform from normalized data into real-valued data. 
#
# @title The data de-normalization
# @param dt.norm a matrix (n x m) of the normalized data.
# @param range.data a matrix (2 x n) containing the range of the data, where n is the number of variables, and
# first and second rows are the minimum and maximum value, respectively. 
# @param min.scale the minimum value within normalization.
# @param max.scale the maximum value within normalization.
# @seealso \code{\link{norm.data}}
# @return the real-valued data
denorm.data <- function(dt.norm, range.data, min.scale = 0, max.scale = 1){
	
	func <- function(i, j, dt.norm, range.data, min.scale, max.scale){
		return(range.data[1, j] + ((dt.norm[i, j] - min.scale) * (range.data[2, j] - range.data[1, j])) / (max.scale - min.scale))
	}
	
	VecFun <- Vectorize(func, vectorize.args=list("i","j"))
	data.denorm <- outer(1 : nrow(dt.norm), 1 : ncol(dt.norm), VecFun, dt.norm, range.data, min.scale, max.scale)
	
	
	return(data.denorm)
}

# This function is to transform from real-valued data into normalized data. 
#
# @title The data normalization
# @param dt.ori a matrix (n x m) of the original data.
# @param range.data a matrix (2 x n) containing the range of the data, where n is the number of variables, and
# first and second rows are the minimum and maximum value, respectively. 
# @param min.scale the minimum value within normalization.
# @param max.scale the maximum value within normalization.
# @seealso \code{\link{denorm.data}}
# @return the normalized data
norm.data <- function(dt.ori, range.data, min.scale = 0, max.scale = 1){

	func <- function(i, j, dt.ori, range.data, min.scale, max.scale){
		return(min.scale + (dt.ori[i, j] - range.data[1, j]) * (max.scale - min.scale) / (range.data[2, j] - range.data[1, j]))
	}
	
	VecFun <- Vectorize(func, vectorize.args=list("i","j"))
	data.norm <- outer(1 : nrow(dt.ori), 1 : ncol(dt.ori), VecFun, dt.ori, range.data, min.scale, max.scale)

	return(data.norm)
}
