#
#     Description of this R script:
#     R routine for creating a sgl algorithm configuration
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

#' Create a new algorithm configuration
#'
#' With the exception of \code{verbose} it is not recommended to change any of the default values.
#' 
#' @param tolerance_penalized_main_equation_loop tolerance threshold.
#' @param tolerance_penalized_inner_loop_alpha tolerance threshold.
#' @param tolerance_penalized_inner_loop_beta tolerance threshold.
#' @param tolerance_penalized_middel_loop_alpha tolerance threshold.
#' @param tolerance_penalized_outer_loop_alpha tolerance threshold.
#' @param tolerance_penalized_outer_loop_beta tolerance threshold.
#' @param tolerance_penalized_outer_loop_gamma tolerance threshold.
#' @param use_bound_optimization if \code{TRUE} hessian bound check will be used.
#' @param use_stepsize_optimization_in_penalizeed_loop if \code{TRUE} step-size optimization will be used.
#' @param stepsize_opt_penalized_initial_t initial step-size.
#' @param stepsize_opt_penalized_a step-size optimization parameter.
#' @param stepsize_opt_penalized_b step-size optimization parameter.
#' @param inner_loop_convergence_limit inner loop convergence limit.
#' @param verbose If \code{TRUE} some information, regarding the status of the algorithm, will be printed in the R terminal.
#' @return A configuration.
#' @examples
#' config.no_progressbar <- sgl.algorithm.config(verbose = FALSE)
#' @author Martin Vincent
#' @export
sgl.algorithm.config <- function(tolerance_penalized_main_equation_loop = 1e-10, 
		tolerance_penalized_inner_loop_alpha = 1e-4, 
		tolerance_penalized_inner_loop_beta = 1, 
		tolerance_penalized_middel_loop_alpha = 0.01, 
		tolerance_penalized_outer_loop_alpha = 0.01, 
		tolerance_penalized_outer_loop_beta = 0, 
		tolerance_penalized_outer_loop_gamma = 1e-5, 
		use_bound_optimization = TRUE, 
		use_stepsize_optimization_in_penalizeed_loop = TRUE, 
		stepsize_opt_penalized_initial_t = 1,
		stepsize_opt_penalized_a = 0.1, 
		stepsize_opt_penalized_b = 0.1,
		inner_loop_convergence_limit = 1e4,
		verbose = TRUE) {
	
	config <- list()
	
	config$tolerance_penalized_main_equation_loop <- tolerance_penalized_main_equation_loop
	
	config$tolerance_penalized_inner_loop_alpha <- tolerance_penalized_inner_loop_alpha
	config$tolerance_penalized_inner_loop_beta <- tolerance_penalized_inner_loop_beta
	
	config$tolerance_penalized_middel_loop_alpha <- tolerance_penalized_middel_loop_alpha
	
	config$tolerance_penalized_outer_loop_alpha <- tolerance_penalized_outer_loop_alpha
	config$tolerance_penalized_outer_loop_beta <- tolerance_penalized_outer_loop_beta	
	config$tolerance_penalized_outer_loop_gamma <- tolerance_penalized_outer_loop_gamma	
	
	config$use_bound_optimization <- use_bound_optimization
	
	config$use_stepsize_optimization_in_penalizeed_loop <- use_stepsize_optimization_in_penalizeed_loop
	config$stepsize_opt_penalized_initial_t <- stepsize_opt_penalized_initial_t
	config$stepsize_opt_penalized_a <- stepsize_opt_penalized_a
	config$stepsize_opt_penalized_b <- stepsize_opt_penalized_b
	
	config$inner_loop_convergence_limit <- as.integer(inner_loop_convergence_limit)
	
	config$verbose <- verbose
	
	return(config)
}

#' Standard algorithm configuration
#'
#' \code{sgl.standard.config <- sgl.algorithm.config()}
#' 
#' @author Martin Vicnet
#' @export
sgl.standard.config <- sgl.algorithm.config();

#' Featch information about the C side configuration of the package
#' @return list
#' 
#' @author Martin Vicnet
#' @useDynLib sglOptim r_pkg_c_config
#' @export
sgl.c.config <- function() {
	.Call("r_pkg_c_config")
}

#' Simulated data set
#'
#' This data set is for testing only.
#'
#' @name test.data
#' @docType data
#' @keywords data
NULL
