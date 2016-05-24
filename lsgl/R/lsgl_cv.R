#
#     Description of this R script:
#     R interface for linear multiple output sparse group lasso routines.
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

#' @title Linear multiple output cross validation using multiple possessors 
#' 
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response matrix, matrix of size \eqn{N \times K}.
#' @param intercept should the model include intercept parameters.
#' @param weights sample weights, vector of size \eqn{N \times K}.
#' @param grouping grouping of features, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the feature. 
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than 1 and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' @param cv.indices a list of indices of a cross validation splitting.
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads the maximal number of threads to be used.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return
#' \item{Yhat}{the cross validation estimated response matrix}
#' \item{Y.true}{the true response matrix, this is equal to the argument \code{y}}
#' \item{cv.indices}{the cross validation splitting used}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models.}
#' @examples
#' 
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 
#' 
#' N <- 100 #number of samples
#' p <- 50 #number of features
#' K <- 25  #number of groups
#' 
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 
#' X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)
#' 
#' ##Do cross validation
#' lambda <- lsgl.lambda(X1, Y1, alpha = 1, d = 15L, lambda.min = 5, intercept = FALSE)
#' fit.cv <- lsgl.cv(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#' 
#' ## Cross validation errors (estimated expected generalization error)
#' Err(fit.cv)
#' @author Martin Vincent
#' @useDynLib lsgl, .registration=TRUE
#' @export
#' @importFrom utils packageVersion
lsgl.cv <- function(x, y, intercept = TRUE,
		weights = NULL,
		grouping = factor(1:ncol(x)), 
		groupWeights = c(sqrt(ncol(y)*table(grouping))),
		parameterWeights =  matrix(1, nrow = ncol(y), ncol = ncol(x)), 
		alpha = 1, lambda, fold = 10L, cv.indices = list(), max.threads = 2L,
		algorithm.config = lsgl.standard.config) 
{
	
	if(!is.matrix(y)) {
		y <- as.matrix(y)
	}
	
	if(nrow(x) != nrow(y)) {
		stop("x and y must have the same number of rows")
	}
	
	if(!is.null(weights)) {
		if(!all(dim(y) == dim(weights))) {
			stop("w and y must have the same dimensions")
		}
	}
	
	# Get call
	cl <- match.call()
	
	# cast
	grouping <- factor(grouping)
	
	# add intercept
	if(intercept) {
		x <- cBind(Intercept = rep(1, nrow(x)), x)
		groupWeights <- c(0, groupWeights)
		parameterWeights <- cbind(rep(0, ncol(y)), parameterWeights)
		grouping <- factor(c("Intercept", as.character(grouping)), levels = c("Intercept", levels(grouping)))
	}
	
	# create data
	group.names <- if(is.null(colnames(y))) 1:ncol(y) else colnames(y)
	data <- create.sgldata(x, y, weights = weights, group.names = group.names)
	
	# Print info
	if(algorithm.config$verbose) {
		
		n_fold <- if(length(cv.indices) == 0) fold else length(cv.indices)
		
		cat("\nRunning lsgl", n_fold, "fold cross validation ")
		if(data$sparseX & data$sparseY) {
			cat("(sparse design and response matrices)")
		}
		if(data$sparseX & !data$sparseY) {
			cat("(sparse design matrix)")
		}
		if(!data$sparseX & data$sparseY) {
			cat("(sparse response matrix)")
		}
		
		cat("\n\n")
		
		print(data.frame('Samples: ' = print_with_metric_prefix(nrow(x)), 
						'Features: ' = print_with_metric_prefix(data$n.covariate), 
						'Models: ' = print_with_metric_prefix(ncol(y)), 
						'Groups: ' = print_with_metric_prefix(length(unique(grouping))), 
						'Parameters: ' = print_with_metric_prefix(length(parameterWeights)),
						check.names = FALSE), 
				row.names = FALSE, digits = 2, right = TRUE)
		cat("\n")
	}
	
	
	# call SglOptimizer function
	if(!is.null(weights)) {
		obj <- "lsgl_w_"
	} else {
		obj <- "lsgl_"
	}
	
	callsym <- paste(obj, if(data$sparseX) "xs_" else "xd_", if(data$sparseY) "ys" else "yd", sep = "")
	
	res <- sgl_cv(callsym, "lsgl", data, grouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, algorithm.config)
	
	# Add weights
	res$weights <- weights
	
	# Add true response
	res$Y.true <- y
	
	# Responses
	res$Yhat <- lapply(res$responses$link, t)
	res$responses <- NULL
	
	res$lsgl_version <- packageVersion("lsgl")
	res$intercept <- intercept
	res$call <- cl
	
	class(res) <- "lsgl"
	return(res)
}
