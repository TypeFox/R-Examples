#
#     Description of this R script:
#     R interface to routine for computing a lambda sequence for the regularization path
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#


#' Generic routine for computing a lambda sequence for the regularization path
#' 
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence.
#' @param lambda.min the smallest lambda value in the computed sequence. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return a vector of length \code{d} containing the compute lambda sequence.
#' @author Martin Vincent
#' @export
#' @import Matrix
sgl_lambda_sequence <- function(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, d, lambda.min, algorithm.config = sgl.standard.config) {
	
	if(lambda.min <= 0) stop("lambda.min should be larger than zero")
	
	# cast
	d <- as.integer(d)
		
	args <- prepare.args(data, parameterGrouping, groupWeights, parameterWeights, alpha)
		
	call_sym <- paste(module_name, "sgl_lambda", sep="_")
	res <- .Call(call_sym, PACKAGE = PACKAGE, args$data, args$block.dim, args$groupWeights, args$parameterWeights, args$alpha, d, lambda.min, algorithm.config)
	
	if(res[1] < res[length(res)]) stop(paste("lamdba.min should be smaller than lambda.max (=", round(res[1],4),")", sep=""))
	
	return(res)
}
