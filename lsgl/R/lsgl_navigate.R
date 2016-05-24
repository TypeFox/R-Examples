#
#     Description of this R script:
#     R interface for linear multi-response models using sparse group lasso.
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

#' @title Compute error rates 
#' 
#' @description
#' Compute and return an error for each model. The error may be spicifed in the \code{loss} argument.
#' 
#' The root-mean-square error (RMSE) is
#' \deqn{\frac{1}{K}\sum_{i = 1}^K \sqrt{\frac{1}{N}\sum_{j=1}^N Y_{ji}-(X\hat \beta)_{ji}}}
#' RMSE is the default error.
#' 
#' The objective value error (OVE) is
#' \deqn{\|Y - X\hat \beta\|_F}
#' 
#' The scaled objective value error (SOVE) is
#' \deqn{\frac{1}{NK}\|Y - X\hat \beta\|_F}
#' 
#' @param object a lsgl object.
#' @param data a design matrix (the \eqn{X} matrix).
#' @param response redirected to \code{y}.
#' @param y a matrix of the true responses (the \eqn{Y} matrix).
#' @param loss the loss (error) function. Either a function taking two arguments or 
#' one of the following character strings \code{RMSE}, \code{OVE} or \code{SOVE}.
#' @param ... ignored.
#' @return a vector of errors.
#' 
#' @author Martin Vincent
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 
#' 
#' N <- 100 #number of samples
#' p <- 50 #number of features
#' K <- 15  #number of groups
#' 
#' # simulate beta matrix and X matrix
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 
#' X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)
#' 
#' X2<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y2 <-X2%*%B+matrix(rnorm(N*K,0,1),N,K)
#' 
#' #### Fit models using X1
#' lambda <- lsgl.lambda(X1, Y1, alpha = 1, d = 25L, lambda.min = 5, intercept = FALSE)
#' fit <- lsgl(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#' 
#' ## Training errors:
#' Err(fit, X1)
#' 
#' ## Errors predicting Y2:
#' Err(fit, X2, Y2)
#' 
#' #### Do cross validation
#' fit.cv <- lsgl.cv(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#' 
#' ## Cross validation errors (estimated expected generalization error)
#' Err(fit.cv)
#'
#' ## Cross validation errors using objective value error measures
#' Err(fit.cv, loss = "OVE")
#'
#' @method Err lsgl
#' @import sglOptim
#' @export
Err.lsgl <- function(object, data = NULL, response = object$Y.true, y = response, loss = "RMSE", ... ) {
	
	if(!is.function(loss)) {
		loss <- switch(loss, 
			RMSE = function(x,y) mean(sqrt(colMeans((x - y)^2))), 
			OVE = function(x,y) sqrt(sum((x - y)^2)),
			SOVE = function(x,y) 1/length(x)*sqrt(sum((x - y)^2)),
			stop("Unknown loss"))
	}

	return(compute_error(object, data = data, response.name = "Yhat", response = y, loss = loss))
}

#' @title Nonzero features
#' 
#' @description
#' Extracts the nonzero features for each model.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return a list of of length \code{nmod(x)} containing the nonzero features (that is nonzero columns of the beta matrices)
#' 
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 
#' 
#' N <- 100 #number of samples
#' p <- 50 #number of covariates
#' K <- 25  #number of groups
#' 
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 
#' 
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	
#' 
#' lambda<-lsgl.lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#' 
#' # the nonzero features of model 1, 10 and 25
#' features(fit)[c(1,10,25)]
#'
#' # count the number of nonzero features in each model
#' sapply(features(fit), length)
#'
#' @author Martin Vincent
#' @method features lsgl
#' @import sglOptim
#' @export
features.lsgl <- function(object, ...) {
	
	class(object) <- "sgl" # Use std function

	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}
	
	return(features(object))
}

#' @title Nonzero parameters
#' 
#' @description
#' Extracts the nonzero parameters for each model.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#' 
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
#' 
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	
#' 
#' lambda<-lsgl.lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#' 
#' # the nonzero parameters of model 1, 10 and 25
#' parameters(fit)[c(1,10,25)]
#'
#' # count the number of nonzero parameters in each model
#' sapply(parameters(fit), sum)
#'
#' @author Martin Vincent
#' @method parameters lsgl
#' @import sglOptim
#' @export
parameters.lsgl <- function(object, ...) {
	
	class(object) <- "sgl" # Use std function
	
	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}
	
	return(parameters(object))
}

#' @title Returns the number of models in a lsgl object
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return the number of models in \code{object}
#' 
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
#' 
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	
#' 
#' lambda<-lsgl.lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#' 
#' # the number of models
#' nmod(fit)
#'
#' @author Martin Vincent
#' @method nmod lsgl
#' @import sglOptim
#' @export
nmod.lsgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(nmod(object, ...))
}

#' @title Exstract the fitted models 
#' 
#' @description
#' Returns the fitted models, that is the estimated \eqn{\beta} matrices.
#' 
#' @param object a lsgl object 
#' @param index indices of the models to be returned
#' @param ... ignored
#' @return a list of \eqn{\beta} matrices.
#' 
#' @author Martin Vincent
#' @method models lsgl
#' @import sglOptim
#' @export
models.lsgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	
	return(models(object, ...))
}

#' @title Extract nonzero coefficients 
#'
#' @param object a lsgl object
#' @param index indices of the models
#' @param ... ignored
#' @return a list of length \code{length(index)} with nonzero coefficients of the models
#'
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 
#' 
#' N <- 100 #number of samples
#' p <- 50 #number of covariates
#' K <- 25  #number of groups
#' 
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 
#' 
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	
#' 
#' lambda<-lsgl.lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # the nonzero coefficients of the models 1, 2 and 5
#' coef(fit, index = c(1,2,5))
#'
#' @author Martin Vincent
#' @method coef lsgl
#' @import sglOptim
#' @importFrom stats coef
#' @export
coef.lsgl <- function(object, index = 1:nmod(object), ...) {
	
	class(object) <- "sgl" # Use std function
	
	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}
	
	return(coef(object, index = index, ...))
}


#' Print function for lsgl
#'
#' This function will print some general information about the lsgl object
#'  
#' @param x lsgl object
#' @param ... ignored
#' 
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 
#' 
#' N <- 100 #number of samples
#' p <- 25 #number of features
#' K <- 15  #number of groups
#' 
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 
#' 
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	
#' 
#' lambda<-lsgl.lambda(X,Y, alpha=1, d = 25, lambda.min= 5, intercept=FALSE)
#' fit <-lsgl(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # Print some information about the estimated models
#' fit
#' 
#' ## Cross validation
#' fit.cv <- lsgl.cv(X, Y, alpha = 1, lambda = lambda, intercept = FALSE)
#' 
#'# Print some information
#' fit.cv
#' 
#' @method print lsgl
#' @author Martin Vincent
#' @import sglOptim
#' @export
print.lsgl <- function(x, ...) {
	sgl_print(x)
}
