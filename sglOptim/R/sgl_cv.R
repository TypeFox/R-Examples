#
#     Description of this R script:
#     R interface generic sparse group lasso cross validation
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

#' Generic sparse group lasso cross validation using multiple possessors 
#' 
#' 
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than \eqn{1} and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' If \code{fold}\eqn{\le}\code{max(table(classes))} then the data will be split into \code{fold} disjoint subsets keeping the ration of classes approximately equal.
#' Otherwise the data will be split into \code{fold} disjoint subsets without keeping the ration fixed.
#' @param cv.indices a list of indices of a cross validation splitting. 
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads the maximal number of threads to be used.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return
#' \item{responses}{content will depend on the C++ response class}
#' \item{cv.indices}{the cross validation splitting used}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models}
#' \item{lambda}{the lambda sequence used.}
#' @author Martin Vincent
#' @importFrom utils packageVersion
#' @export
sgl_cv <- function(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, lambda, fold = 2, cv.indices = list(), max.threads = 2, algorithm.config = sgl.standard.config) {
	
	# cast
	fold <- as.integer(fold)
		
	if(length(cv.indices) == 0) {
		
		# Check fold
		if(fold < 2) {
			stop("fold must be equal to or larger than 2")
		}
		
		if(fold > length(data$G)) {
			stop("fold must be equal to or less than the number of samples")
		}
		
		if(fold > min(table(data$G)) || length(unique(data$G)) == 1) {
			# use random sample indices
			use.cv.indices <- TRUE
			cv.indices <- split(sample(1:(length(data$G))), 1:fold)
			# TODO need to ensure that each split has one sample from each class
			
		} else {
			# compute cv indices
			cv.indices <- lapply(unique(data$G), function(x) .divide_indices(which(data$G == x), fold))
			cv.indices <- lapply(cv.indices, function(x) sample(x))
			cv.indices <- lapply(1:fold, function(i) sort(unlist(lapply(cv.indices, function(x) x[[i]]))))
		}
		
		# Chek consistency of cv.indices
		if(length(unlist(cv.indices)) != length(data$G) || sum(duplicated(unlist(cv.indices))) > 0) {
			stop("Internal error computing the cross validation splitting (this is a bug)")
		}
		
	} else {
		# use user supplied cv splitting
		# Chek consistency of cv.indices
		if(length(unlist(cv.indices)) != length(data$G) || sum(duplicated(unlist(cv.indices))) > 0) {
			stop("User supplied cv.indices are invalid (the cv.indices does not represent a cross validation splitting, use subsampling for general subsampling)")
		}
	}

	samples <- 1:max(unlist(cv.indices)) 
	training <- lapply(cv.indices, function(x) samples[-x])
	test <- cv.indices
	
	res <- sgl_subsampling(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, TRUE, max.threads, algorithm.config)
	
	# Add cv.indices
	res$cv.indices <- cv.indices
	
	# Set version, type and class and return
	res$sglOptim_version <- packageVersion("sglOptim")
	res$type <- "cv"
	class(res) <- "sgl"
	return(res)
}

.divide_indices <- function(indices, fold) {
	
	if(fold > length(indices)) {
		stop("fold large than length of indices vector")
	}
	
	tmp <- sapply(0:(length(indices)-1), function(i) (i %% fold)+1)	
	return(lapply(1:fold, function(i) indices[tmp == i]))
}