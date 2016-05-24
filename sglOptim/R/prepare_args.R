#
#     Description of this R script:
#     Routines for handling sgl-data
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


#' @title Generic rearrange function
#'
#' @description
#' Rearrange the order of the covariates in the \code{data} object.
#'
#' @param data a data object
#' @param covariate.order the new order of the covariates
#' @param ... additional parameters
#' @return a rearranged data object of same class as \code{data}
#' 
#' @seealso rearrange.sgldata
#' @author Martin Vincent
#' @export
rearrange <- function(data, covariate.order, ...) UseMethod("rearrange")

#' @title Generic function for preparing the sgl call arguments
#' 
#' @description
#' Compute and prepare the sgl call arguments for the objective function
#' \deqn{\mathrm{loss}(\mathrm{data})(\beta) + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\mathrm{loss}} is a loss/objective function.
#' The \eqn{n} parameters are organized in the parameter matrix \eqn{\beta} with dimension \eqn{q\times p}.
#' The vector \eqn{\beta^{(J)}} denotes the \eqn{J} parameter group, the dimension of \eqn{\beta^{(J)}} is denote by \eqn{d_J}.
#' The dimensions \eqn{d_J} must be multiple of \eqn{q}, and \eqn{\beta = (\beta^{(1)} \cdots \beta^{(m)})}.
#' The group weights \eqn{\gamma \in [0,\infty)^m} and the parameter weights \eqn{\xi \in [0,\infty)^{qp}}.
#' @param data a data object
#' @param ... additional parameters
#' @return 
#' \item{block.dim}{a vector of length \eqn{m}, containing the dimensions \eqn{d_J} of the groups (i.e. the number of parameters in the groups)}
#' \item{groupWeights}{a vector of length \eqn{m}, containing the group weights}
#' \item{parameterWeights}{a matrix of dimension \eqn{q \times p}, containing the parameter weights}
#' \item{alpha}{the \eqn{\alpha} value}
#' \item{data}{the data parsed to the loss module}
#' \item{group.order}{original order of the columns of \eqn{\beta}. Before sgl routines return \eqn{\beta} will be reorganized according to this order.}
#'
#' @seealso prepare.args.sgldata
#' @author Martin Vincent
#' @export
#' @family sgldata
prepare.args <- function(data, ...) UseMethod("prepare.args")
	

#' @title Rearrange sgldata
#'
#' @description
#' Rearrange the order of the covariates in a sgldata object.
#'
#' @param data  a sgldata object
#' @param covariate.order the new order of the covariates
#' @param ... not used
#' @return a sgldata object with the covariates reordered
#' @author Martin Vincent
#' @method rearrange sgldata
#' @export
#' @family sgldata
rearrange.sgldata <- function(data, covariate.order, ...) 
{
	
	data$X <- data$X[,covariate.order, drop = FALSE]
	data$covariate.names <- data$covariate.names[covariate.order]
	
	return(data)
}

#' @title Create a sgldata object
#'
#' @description
#' Creates a sgldata object from a design matrix and an optional response vector or matrix.
#'
#' @param x the design matrix, a matrix of size \eqn{N \times p} (will be parsed to the loss module as \code{X}).
#' @param y the responses, \code{NULL}, a vector or a matrix (will be parsed to the loss module as \code{Y})..
#' @param weights sample weights, a vector of length \eqn{N} (will be parsed to the loss module as \code{W}).
#' @param sampleGrouping grouping of samples, a factor of length \eqn{N} (will be parsed to the loss module as \code{G}). Default is no grouping (NULL), that is all samples is the same group.
#' @param group.names a vector with the names of the parameter groups (the length must equal the number of rows in the \eqn{\beta} matrix). 
#' @param sparseX if TRUE \code{x} will be treated as sparse, if FALSE \code{x} will be treated as dens.
#' @param sparseY if TRUE \code{y} will be treated as sparse, if FALSE \code{y} will be treated as dens.
#' @author Martin Vincent
#' @importFrom methods is
#' @importFrom methods as
#' @export
#' @family sgldata
create.sgldata <- function(x, y, weights = NULL, sampleGrouping = NULL, group.names = NULL, sparseX = is(x, "sparseMatrix"), sparseY = is(y, "sparseMatrix")) {
	
	#TODO dim checks
	
	### Data 
	data <- list()
	
	### Dimensions 
	data$n.covariate <- ncol(x)
	data$n.samples <- nrow(x)		
	
	### X data
	data$sparseX <- sparseX

	if(data$sparseX) {
		data$X <- as(x, "CsparseMatrix")
	} else if(is(x, "matrix") || is(x, "sparseMatrix")){
		data$X <- as.matrix(x)
	} else {
		stop("design matrix of unkown type")
	}

	# Dim-names
	data$sample.names <- rownames(x)
	data$covariate.names <- colnames(x)
	
	### Y data
	data$sparseY <- sparseY
		
	if(is.null(y)) {
		data$Y <- NULL
	} else if(is.vector(y)) {
		data$Y <- as.numeric(y)
	} else if(is.matrix(y)) {
		data$Y <- apply(y, 2, as.numeric)
	} else if(sparseY) {
		data$Y <- as(y, "CsparseMatrix")
	} else {
		stop("y is of unknown type")
	}
		
	### weights
	if(is.null(weights)) {
		data$W <- rep(1/data$n.samples, data$n.samples)
	} else if(is.vector(weights)) {
		data$W <- as.numeric(weights)
	} else if(is.matrix(weights)) {
		data$W <- apply(weights, 2, as.numeric)
	} else {
		stop("weights is of unknown type")
	}

	### sample grouping
	if(is.null(sampleGrouping)) {
		sampleGrouping <- rep(1, data$n.samples)
	}
	
	sampleGrouping <- factor(sampleGrouping)
	data$G <- as.integer(factor(sampleGrouping))-1L
	
	### group names 
	if(is.null(group.names)) {
		data$group.names <- levels(sampleGrouping)
	} else {
		data$group.names <- group.names
	}
	
	data$n.groups <- as.integer(length(data$group.names))
	
	class(data) <- "sgldata"
	return(data)
}

#' @title Prepare sgl function arguments 
#' 
#' @description
#' Prepare sgl function arguments using sgldata.
#'
#' @param data a sgldata object
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}, that is the same dimension as \eqn{\beta}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param ... not used
#' @method prepare.args sgldata
#' @export
#' @family sgldata
#' @author Martin Vincent
prepare.args.sgldata <- function(data, parameterGrouping, groupWeights, parameterWeights, alpha, ...) {
	
	# If Lasso then ignore grouping
	if(alpha == 1) {
		parameterGrouping <- factor(1:data$n.covariate)
		groupWeights <- rep(1, data$n.covariate)
	}
	
	ncov <- data$n.covarites
	ngrp <- data$n.groups
		
	group.order <- order(parameterGrouping)
	
	### Reorder data
	data <- rearrange(data, group.order)
	
	parameterWeights <- parameterWeights[,group.order, drop = FALSE]
			
	### Compute block dim
	block.dim <- ngrp*as.integer(table(parameterGrouping))
	
	### Prepare data format
	# sparse X format
	if(data$sparseX) {
		data$X <- list(dim(data$X), data$X@p, data$X@i, data$X@x)
	}
	
	# sparse Y format
	if(data$sparseY) {
		data$Y <- list(dim(data$Y), data$Y@p, data$Y@i, data$Y@x)
	}
	
	### Create args list
	args <- list()
	args$block.dim <- block.dim
	args$groupWeights <- groupWeights
	args$parameterWeights <- parameterWeights
	args$alpha <- alpha
	args$data <- data
	args$group.order <- group.order
	
	return(args)
}
