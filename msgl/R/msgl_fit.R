#
#     Description of this R script:
#     R interface for multinomial sparse group lasso rutines.
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

#' @title Fit a multinomial sparse group lasso regularization path.
#'
#' @description 
#' Fit a sequence of multinomial logistic regression models using sparse group lasso, group lasso or lasso.
#' In addition to the standard parameter grouping the algorithm supports further grouping of the features.
#'
#' @details 
#' For a classification problem with  \eqn{K} classes and \eqn{p} features (covariates) dived into \eqn{m} groups.
#' This function computes a sequence of minimizers (one for each lambda given in the \code{lambda} argument) of
#' \deqn{\hat R(\beta) + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\hat R} is the weighted empirical log-likelihood risk of the multinomial regression model.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of features
#' (default is one covariate per group, hence the default dimension of \eqn{\beta^{(J)}} is \eqn{K}).
#' The group weights \eqn{\gamma \in [0,\infty)^m} and parameter weights \eqn{\xi \in [0,\infty)^n} may be explicitly specified.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of features, a vector of length \eqn{p}. Each element of the vector specifying the group of the feature.
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups).
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and 
#' \deqn{\sqrt{K\cdot\textrm{number of features in the group}}}
#' for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times p}.
#' If \code{parameterWeights = NULL} default weights will be used.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the features are standardize before fitting the model. The model parameters are returned in the original scale.
#' @param lambda the lambda sequence for the regularization path.
#' @param return the indices of lambda values for which to return a the fitted parameters.
#' @param intercept should the model fit include intercept parameters (note that due to standardization the returned beta matrix will always have an intercept column)
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param algorithm.config the algorithm configuration to be used.
#' @return
#' \item{beta}{the fitted parameters -- a list of length \code{length(lambda)} with each entry a matrix of size \eqn{K\times (p+1)} holding the fitted parameters}
#' \item{loss}{the values of the loss function}
#' \item{objective}{the values of the objective function (i.e. loss + penalty)}
#' \item{lambda}{the lambda values used}
#' \item{classes.true}{the true classes used for estimation, this is equal to the \code{classes} argument}
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Model 10, i.e. the model corresponding to lambda[10]
#' models(fit)[[10]] 
#' 
#' # The nonzero features of model 10
#' features(fit)[[10]]
#' 
#' # The nonzero parameters of model 10
#' parameters(fit)[[10]]
#' 
#' # The training errors of the models.
#' Err(fit, x) 
#' # Note: For high dimensional models the training errors are almost always over optimistic,
#' # instead use msgl.cv to estimate the expected errors by cross validation 
#' 
#' @author Martin Vincent
#' @importFrom utils packageVersion
#' @importFrom methods is
#' @export
#' @useDynLib msgl, .registration=TRUE
#' @import Matrix
msgl <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE,
		lambda, return = 1:length(lambda), intercept = TRUE, sparse.data = is(x, "sparseMatrix"), algorithm.config = msgl.standard.config) {
	
	# Get call
	cl <- match.call()
	
	#Check dimensions 
	if(nrow(x) != length(classes)) {
		stop("the number of rows in x must match the length of classes")
	}
	
	# Default values
	if(is.null(grouping)) {
		covariateGrouping <- factor(1:ncol(x))
	} else {
		# ensure factor
		covariateGrouping <- factor(grouping)
	}
	
	# cast
	classes <- factor(classes)
	return <- as.integer(return)
	
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(length(levels(classes))*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = length(levels(classes)), ncol = ncol(x))
	}
	
	# Standardize
	if(standardize) {
		
		if(sparse.data) {
			x.scale <- sqrt(colMeans(x*x) - colMeans(x)^2)
			x.center <- rep(0, length(x.scale))
			x <- x%*%Diagonal(x=1/x.scale)
		} else {
			x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
			x.scale <- attr(x, "scaled:scale")
			x.center <- if(sparse.data) rep(0, length(x.scale)) else attr(x, "scaled:center")
		}
	}
	
	if(intercept) {
		intercept.value = 1
	} else {
		intercept.value = 0
	}
		# add intercept
	if(is.null(colnames(x))) {
		x <- cBind(rep(intercept.value, nrow(x)), x)
	} else {
		x <- cBind(Intercept = rep(intercept.value, nrow(x)), x)
	}
	
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(classes))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))

	# create data
	data <- create.sgldata(x, y = NULL, sampleWeights, classes, sparseX = sparse.data)
	
	# call sglOptim function
	if(data$sparseX) {
		
		if(algorithm.config$verbose) {
			
			cat("\nRunning msgl (sparse design matrix)\n\n")
			print(data.frame('Samples: ' = print_with_metric_prefix(length(sampleWeights)), 
							'Features: ' = print_with_metric_prefix(data$n.covariate), 
							'Classes: ' = print_with_metric_prefix(length(levels(classes))), 
							'Groups: ' = print_with_metric_prefix(length(unique(covariateGrouping))), 
							'Parameters: ' = print_with_metric_prefix(length(parameterWeights)),
							check.names = FALSE), 
					row.names = FALSE, digits = 2, right = TRUE)
			cat("\n")
		}
		
		res <- sgl_fit("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)
	} else {
		
		if(algorithm.config$verbose) {
			cat("\nRunning msgl (dense design matrix) \n\n")
			print(data.frame('Samples: ' = print_with_metric_prefix(length(sampleWeights)), 
							'Features: ' = print_with_metric_prefix(data$n.covariate), 
							'Classes: ' = print_with_metric_prefix(length(levels(classes))), 
							'Groups: ' = print_with_metric_prefix(length(unique(covariateGrouping))), 
							'Parameters: ' = print_with_metric_prefix(length(parameterWeights)),
							check.names = FALSE), 
					row.names = FALSE, digits = 2, right = TRUE)
			cat("\n")
		}
		
		res <- sgl_fit("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)
	}
	
	# Convert beta back to the org scale
	if(standardize) {
		res$beta <- .to_org_scale(beta = res$beta, x.scale = x.scale, x.center = x.center)
	}
	
	res$classes.true <- classes

	# Various 
	res$msgl_version <- packageVersion("msgl")
	res$call <- cl

	class(res) <- "msgl"
	return(res)
}

#' @title Computes a lambda sequence for the regularization path
#'
#' @description
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of features, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate.
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups).
#' The first element of the vector is the intercept weight.
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of features in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}.
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence
#' @param standardize if TRUE the features are standardize before fitting the model. The model parameters are returned in the original scale.
#' @param lambda.min the smallest lambda value in the computed sequence.
#' @param intercept should the model include intercept parameters
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param algorithm.config the algorithm configuration to be used.
#' @return a vector of length \code{d} containing the computed lambda sequence.
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#'
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100, lambda.min = 0.01)
#' @author Martin Vincent
#' @importFrom methods is
#' @export
#' @useDynLib msgl, .registration=TRUE
msgl.lambda.seq <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, d = 100L, standardize = TRUE, lambda.min, intercept = TRUE, sparse.data = is(x, "sparseMatrix"), algorithm.config = sgl.standard.config) {

		#Check dimensions 
		if(nrow(x) != length(classes)) {
			stop("the number of rows in x must match the length of classes")
		}
	
        # cast
        classes <- factor(classes)

        # Default values
        if(is.null(grouping)) {
                covariateGrouping <- factor(1:ncol(x))
        } else {
                # ensure factor
                covariateGrouping <- factor(grouping)
        }

        if(is.null(groupWeights)) {
                groupWeights <- c(sqrt(length(levels(classes))*table(covariateGrouping)))
        }

        if(is.null(parameterWeights)) {
                parameterWeights <-  matrix(1, nrow = length(levels(classes)), ncol = ncol(x))
        }

	    # Standardize
        if(standardize) {
			
			if(sparse.data) {
				x.scale <- sqrt(colMeans(x*x) - colMeans(x)^2)
				x.center <- rep(0, length(x.scale))
				x <- x%*%Diagonal(x=1/x.scale)
			} else {
				x <- scale(x, TRUE, TRUE)
                x.scale <- attr(x, "scaled:scale")
                x.center <- attr(x, "scaled:center")
			}
        }

		if(intercept) {
			intercept.value = 1
		} else {
			intercept.value = 0
		}
		# add intercept
		if(is.null(colnames(x))) {
			x <- cBind(rep(intercept.value, nrow(x)), x)
		} else {
			x <- cBind(Intercept = rep(intercept.value, nrow(x)), x)
		}
		
		groupWeights <- c(0, groupWeights)
		parameterWeights <- cbind(rep(0, length(levels(classes))), parameterWeights)
		covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
		
        # create data
        data <- create.sgldata(x, y = NULL, sampleWeights, classes, sparseX = sparse.data)

        # call SglOptimizer function
        if(data$sparseX) {
                lambda <- sgl_lambda_sequence("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
        } else {
                lambda <- sgl_lambda_sequence("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
        }

        return(lambda)
}


.to_org_scale <- function(beta, x.scale, x.center) {

		for(l in 1:length(beta)) {

                beta.org <- t(t(beta[[l]])*c(1,1/x.scale))
                beta.org[,1] <- beta.org[,1] - rowSums(t(t(beta[[l]][,-1])*(x.center/x.scale)))

                beta[[l]] <- beta.org
        }

        return(beta)
}


#' Simulated data set
#'
#' The use of this data set is only intended for testing and examples.
#' The data set contains 100 simulated samples grouped into 10 classes.
#' For each sample 400 features have been simulated.
#'
#' @name sim.data
#' @docType data
#' @keywords data
NULL

