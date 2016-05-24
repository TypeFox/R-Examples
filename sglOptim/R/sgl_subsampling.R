#
#     Description of this R script:
#     R interface to generic sparse group lasso subsampling procedure
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

#' Generic sparse group lasso subsampling procedure
#'
#' Support the use of multiple processors.
#' 
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param training a list of training samples, each item of the list corresponding to a subsample.
#' Each item in the list must be a vector with the indices of the training samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{test} list.  
#' @param test a list of test samples, each item of the list corresponding to a subsample.
#' Each item in the list must be vector with the indices of the test samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{training} list.
#' @param collapse if \code{TRUE} the results for each subsample will be collapse into one result (this is useful if the subsamples are not overlapping)
#' @param max.threads the maximal number of threads to be used.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return
#' \item{responses}{content will depend on the C++ response class}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models}
#' \item{lambda}{the lambda sequence used.}
#' @author Martin Vincent
#' @useDynLib sglOptim, .registration=TRUE
#' @importFrom utils packageVersion
#' @export
sgl_subsampling <- function(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, collapse = FALSE, max.threads = 2, algorithm.config = sgl.standard.config) {
	
	# cast
	max.threads <- as.integer(max.threads)
		
	#Check training and test consistency:
	if(!is.list(training) | !is.list(test)) {
		stop("The arguments traning and test should be lists")
	}
	
	if(length(training) != length(test)) {
		stop("The length of the lists traning and test should be equal")
	}
	
	# Prapare arguments
	args <- prepare.args(data, parameterGrouping, groupWeights, parameterWeights, alpha)
	
	training <- lapply(training, sort)
	test <- lapply(test, sort)
	
	training.0 <- lapply(training, function(x) as.integer(x - 1))
	test.0 <- lapply(test, function(x) as.integer(x - 1))
	
	call_sym <- paste(module_name, "sgl_subsampling", sep="_")
	res <- .Call(call_sym, PACKAGE = PACKAGE, args$data, args$block.dim, args$groupWeights, args$parameterWeights, args$alpha, lambda, training.0, test.0, collapse, max.threads, algorithm.config)
	
	# Sample names
	sample.names <- data$sample.names
	
	if(collapse == TRUE) {
		# Reorder responses and set sample names
		sample.order <- order(unlist(test))
		res$responses <- lapply(res$responses, function(x) .order_response(x, sample.order))
		res$responses <- lapply(res$responses, function(x) .set_sample_names(x, sample.names))
	} else {
		#Set sample names
		res$responses <- lapply(res$responses, function(x) lapply(1:length(x), function(i) .set_sample_names(x[[i]], sample.names[test[[i]]])))
	}
	
	# Names
	rownames(res$features) <- paste("subsample", 1:length(training))
	rownames(res$parameters) <- paste("subsample", 1:length(training))
	
	res$lambda <- lambda
	
	# Set version, type and class and return
	res$sglOptim_version <- packageVersion("sglOptim")
	res$type <- "subsampling"
	class(res) <- "sgl"
	
	return(res)
}

.set_sample_names <- function(response, sample.names) {
	
	if(is.list(response)) {
		return(lapply(response, function(x) t(.set_sample_names(t(x), sample.names))))
	}
	
	if(is.matrix(response)) {
		rownames(response) <- sample.names
		return(response)
	}
	
	stop("Unknown response type")
	
}

.order_response <- function(response, sample.order) {
	
	if(is.list(response)) {
		return(lapply(response, function(x) t(.order_response(t(x), sample.order))))
	}
	
	if(is.matrix(response)) {
		return(response[sample.order,])
	}
	
	stop("Unknown response type")
}
