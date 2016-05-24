#
#     Description of this R script:
#     R interface for multinomial sparse group lasso routines.
#
#     Intended for use with R.
#     Copyright (C) 2013 Martin Vincent
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
#' Compute error rates. 
#' If \code{type = "rate"} then the misclassification rates will be computed.
#' If \code{type = "count"} then the misclassification counts will be computed.
#' If \code{type = "loglike"} then the negative log likelihood error will be computed.
#'
#' @param object a msgl object
#' @param data a matrix of 
#' @param response a vector of classes
#' @param classes a vector of classes
#' @param type type of error rate
#' @param ... ignored
#' @return a vector of error rates
#' 
#' @author Martin Vincent
#' @examples
#' data(SimData)
#' x.all <- sim.data$x
#' x.1 <- sim.data$x[1:50,]
#' x.2 <- sim.data$x[51:100,]
#' classes.all <- sim.data$classes
#' classes.1 <- sim.data$classes[1:50]
#' classes.2 <- sim.data$classes[51:100]
#' 
#' #### Fit models using x.1
#' lambda <- msgl.lambda.seq(x.1, classes.1, alpha = .5, d = 25, lambda.min = 0.075)
#' fit <- msgl(x.1, classes.1, alpha = .5, lambda = lambda)
#' 
#' #### Training errors:
#' 
#' # Misclassification rate
#' Err(fit, x.1)
#' 
#' # Misclassification count
#' Err(fit, x.1, type = "count")
#' 
#' # Negative log likelihood error
#' Err(fit, x.1, type="loglike")
#'  
#' # Misclassification rate of x.2
#' Err(fit, x.2, classes.2)
#' 
#' #### Do cross validation
#' fit.cv <- msgl.cv(x.all, classes.all, alpha = .5, lambda = lambda)
#' 
#' #### Cross validation errors (estimated expected generalization error)
#' 
#' # Misclassification rate
#' Err(fit.cv)
#' 
#' # Negative log likelihood error
#' Err(fit.cv, type="loglike")
#' 
#' #### Do subsampling
#' test <- list(1:20, 21:40)
#' train <- lapply(test, function(s) (1:length(classes.all))[-s])
#'
#' fit.sub <- msgl.subsampling(x.all, classes.all, alpha = .5, 
#'  lambda = lambda, training = train, test = test)
#' 
#' # Mean misclassification error of the tests
#' Err(fit.sub)
#' 
#' # Negative log likelihood error
#' Err(fit.sub, type="loglike")
#'  
#' @method Err msgl
#' @export
#' @import sglOptim
Err.msgl <- function(object, data = NULL, response = object$classes.true, classes = response, type = "rate", ... ) {
	
	if(type=="rate") {
		return(compute_error(object, data = data, response.name = "classes", response = classes, loss = function(x,y) mean(x != y)))
	}
	
	if(type=="count") {
		return(compute_error(object, data = data, response.name = "classes", response = classes, loss = function(x,y) sum(x != y)))
	}
	
	if(type=="loglike") {
		loss <- function(x,y) -mean(log(sapply(1:length(y), function(i) x[as.integer(y[i]),i])))
		return(compute_error(object, data = data, response.name = "response", response = classes, loss = loss))
	}
	
	stop("Unknown type")
	
}

#' @title Nonzero features
#' 
#' @description
#' Extracts the nonzero features for each model.
#'
#' @param object a msgl object
#' @param ... ignored
#' @return a list of of length \code{nmod(x)} containing the nonzero features (that is nonzero colums of the beta matrices)
#' 
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#'
#' # the nonzero features of model 1, 10 and 25
#' features(fit)[c(1,10,25)]
#'
#' # count the number of nonzero features in each model
#' sapply(features(fit), length)
#'
#' @author Martin Vincent
#' @method features msgl
#' @import sglOptim
#' @export
features.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(features(object))
}

#' @title Nonzero parameters
#' 
#' @description
#' Extracts the nonzero parameters for each model.
#'
#' @param object a msgl object
#' @param ... ignored
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#' 
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#'
#' # the nonzero parameters of model 1, 10 and 25
#' parameters(fit)[c(1,10,25)]
#'
#' # count the number of nonzero parameters in each model
#' sapply(parameters(fit), sum)
#'
#' @author Martin Vincent
#' @method parameters msgl
#' @import sglOptim
#' @export
parameters.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(parameters(object))
}

#' @title Returns the number of models in a msgl object
#'
#' @param object a msgl object
#' @param ... not used
#' @return the number of models in \code{object}
#' 
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#'
#' # the number of models
#' nmod(fit)
#'
#' @author Martin Vincent
#' @method nmod msgl
#' @import sglOptim
#' @export
nmod.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(nmod(object, ...))
}

#' @title Exstract the fitted models 
#' 
#' @description
#' Returns the fitted models, that is the estimated \eqn{\beta} matrices.
#' 
#' @param object a msgl object 
#' @param index indices of the models to be returned
#' @param ... ignored
#' @return a list of \eqn{\beta} matrices.
#' 
#' @author Martin Vincent
#' @method models msgl
#' @import sglOptim
#' @export
models.msgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	return(models(object, ...))
}

#' @title Extract nonzero coefficients 
#'
#' @param object a msgl object
#' @param index indices of the models
#' @param ... ignored
#' @return a list of length \code{length(index)} with nonzero coefficients of the models
#'
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#'
#' # the nonzero coefficients of the models 1, 10 and 20
#' coef(fit, index = c(1,10,20))
#'
#' @author Martin Vincent
#' @importFrom stats coef
#' @method coef msgl
#' @import sglOptim
#' @export
coef.msgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	return(coef(object, index = index, ...))
}


#' Print function for msgl
#'
#' This funtion will print some general information about the msgl object
#'  
#' @param x msgl object
#' @param ... ignored
#' 
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' 
#' ### Estimation
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25, lambda.min = 0.075)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#'
#' # Print some information about the estimated models
#' fit
#'
#' ### Cross validation
#' fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Print some information
#' fit.cv
#' 
#' ### Subsampling
#' test <- list(1:20, 21:40)
#' train <- lapply(test, function(s) (1:length(classes))[-s])
#'
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, training = train, test = test)
#' 
#' # Print some information
#' fit.sub
#' 
#' @method print msgl
#' @author Martin Vincent
#' @import sglOptim
#' @export
print.msgl <- function(x, ...) {
	sgl_print(x)
}
