#
#     Description of this R script:
#     Routines for navigating sgl objects
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

# S3 functions:

#' @title Generic function for computing error rates
#' 
#' @description
#' Compute and returns an error rate for each model contained in \code{x}.
#' See details for generic use cases.
#'
#' @details
#' The following generic use case should be supported (see for example \pkg{msgl} package for an implementation):
#' 
#' \enumerate{
#' \item With \code{fit} a sgl fit object with models estimated using \code{x} data, the code 
#' 
#' \code{Err(fit, x)}
#' 
#' should return a vector with the \emph{training errors} of the models.
#' 
#' \item With \code{x.new} a new data set with known responses \code{response.new}, the code 
#' 
#' \code{Err(fit, x.new, response.new)}
#' 
#' should return a vector with the errors of the models when applied to the new data set.
#
#' \item With \code{fit.cv} a sgl cross validation object, the code 
#' 
#' \code{Err(fit.cv)}
#' 
#' should return a vector with estimates of the \emph{expected generalization errors} of the models (i.e. the cross validation errors). 
#' 
#' \item If subsampling is supported then, with \code{fit.sub} a sgl subsampling object, the code  
#' 
#' \code{Err(fit.sub)}
#' 
#' should return a matrix with the test errors (each column corresponding to a model, i.e. rows corresponds to tests).
#' }
#'
#' @seealso compute_error
#' @param object an object
#' @param data a data object 
#' @param response a response object
#' @param ... additional parameters (optional)
#' @return 
#' a vector of length \code{nmod(object)} or a matrix with \code{nmod(object)} columns containing error rates for the models
#' @author Martin Vincent
#' @export
Err <- function(object, data, response, ... ) UseMethod("Err")
	

#' @title Helper function for computing error rates
#' @details
#' This function can be used to compute error rates. 
#' It is consist with the use cases of the \code{Err} genetic function. 
#' (see \pkg{msgl} package for an example of how to use this function)
#' @param object a object containing responses
#' @param data a data object
#' @param response.name the name of the response, if \code{response.name = NULL} then \code{object} will be treated as the response.
#' @param response the response
#' @param loss the loss function
#' @return a vector with the computed error rates
#' 
#' @author Martin Vincent
#' @importFrom stats predict
#' @export
compute_error <- function(object, data = NULL, response.name, response, loss) {
	
	if(!is.null(data)) {
		return(compute_error(object = predict(object, data), data = NULL, response.name = response.name, response = response, loss = loss))
	}
	
	if(is.null(response.name) || any(names(object) == response.name)) {
		
		if(is.null(response.name)) {
			r <- object
		} else {
			r <- object[[response.name]]
		}
		
		if(is.list(r) && is.list(response)) {
			return(t(sapply(1:length(r), function(i) compute_error(object = r[[i]], data = NULL, response.name = NULL, response = response[[i]], loss = loss))))
			
		} else if(is.list(r)) {
			return(sapply(r, function(m) loss(m, response)))	
			
		} else if(is.matrix(r)) {
			return(apply(r, 2, FUN = function(v) loss(v, response)))		
			
		} else if(is.vector(r))  {
			return(loss(r, response))
			
		} else {
			stop("Unknown response type")
		}
	}		

	stop(paste("response '", response.name, "' not found", sep =""))

}

#' @title Generic function for extracting nonzero features (or groups) 
#' 
#' @description
#' Extracts nonzero features for each model.
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return a list of length \code{nmod(x)} containing the nonzero features of the models.
#' 
#' @author Martin Vincent
#' @export
features <- function(object, ...) UseMethod("features")

#' @title Generic function for extracting nonzero parameters
#' 
#' @description
#' Extracts nonzero parameters for each model.
#'
#' @param object an object
#' @param ... additional paramters (optional)
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#' 
#' @author Martin Vincent
#' @export
parameters <- function(object, ...) UseMethod("parameters")


#' @title Generic function for counting the number of models
#' 
#' @description
#' Returns the number of models
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return the number of models contained in the object \code{x}.
#' 
#' @author Martin Vincent
#' @export
nmod <- function(object, ...) UseMethod("nmod")

#' @title Generic function for extracting the fitted models 
#' 
#' @description
#' Returns the fitted models
#'
#' @param object an object
#' @param index a vector of indices of the models to be returned
#' @param ... additional parameters (optional)
#' @return a list of length \code{length(index)} containing the models
#' 
#' @author Martin Vincent
#' @export
models <- function(object, index, ...) UseMethod("models")


#' @title Extracting nonzero features
#' 
#' @param object a sgl object 
#' @param ... not used
#' @return a list of vectors containing the nonzero features (that is nonzero columns of the \eqn{beta} matrices)
#' 
#' @author Martin Vincent
#' @method features sgl
#' @export
features.sgl <- function(object, ...) {
	
	if(is.null(object$beta)) {
		stop("object contains no models")
	}
	
	if(is.null(colnames(object$beta[[1]])) || any(duplicated(colnames(object$beta[[1]])))) {
		res <- lapply(object$beta, function(beta) which(colSums(beta != 0) != 0))
	} else {
		res <- lapply(object$beta, function(beta) colnames(beta)[colSums(beta != 0) != 0])
	}
	
	return(res)
}

#' @title Extracting nonzero parameters
#'
#' @param object a sgl object
#' @param ... not used
#' @return a list of vectors containing the nonzero parameters (that is nonzero entries of the \eqn{beta} matrices)
#' 
#' @author Martin Vincent
#' @method parameters sgl
#' @export
parameters.sgl <- function(object, ...) {
	
	if(is.null(object$beta)) {
		stop("object contains no models")
	}
	
	tmp <- features(object)
	res <- sapply(1:length(object$beta), function(i) object$beta[[i]][,tmp[[i]], drop = FALSE] != 0)
		
	return(res)
}

#' @title Returns the number of models in a sgl object
#'
#' @param object a sgl object
#' @param ... not used
#' @return the number of models in \code{object}
#' 
#' @author Martin Vincent
#' @method nmod sgl
#' @export
nmod.sgl <- function(object, ...) {
	return(length(object$lambda))
}

#' @title Returns the estimated models (that is the \eqn{beta} matrices)
#' 
#' @param object a sgl object
#' @param index indices of the models to be returned
#' @param ... not used
#' @return a list of sparse matrices
#' 
#' @author Martin Vincent
#' @method models sgl
#' @export
models.sgl <- function(object, index = 1:nmod(object), ...) {
	
	if(is.null(object$beta)) {
		stop("object contains no models")
	}
	
	return(object$beta[index])
}

#' @title Extracting the nonzero coefficients 
#'
#' @param object a sgl object
#' @param index indices of the models
#' @param ... not used
#' @return a list of with nonzero coefficients of the models
#' 
#' @author Martin Vincent
#' @method coef sgl
#' @export
coef.sgl <- function(object, index = 1:nmod(object), ...) {
	
	if(is.null(object$beta)) {
		stop("object contains no models")
	}
	
	return(lapply(object$beta[index], function(beta) beta[,colSums(beta != 0) != 0, drop = FALSE]))
}


#' @title Print information about sgl object
#' 
#' @param x a object of sgl family class 
#'  
#' @author Martin Vincent
#' @export
sgl_print <- function(x) {
	
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
			"\n", sep = "")
	
	if(x$type == "fit") {
		
		cat("\nModels:\n\n")
		sel <- 1:5*floor(nmod(x)/5)
		
		feat <- sapply(features(x), length)
		para <- sapply(parameters(x), sum)
		
		print(data.frame('Index: ' = sel, 
						'Lambda: ' = x$lambda[sel], 
						'Features: ' = print_with_metric_prefix(feat[sel]), 
						'Parameters: ' = print_with_metric_prefix(para[sel]), check.names = FALSE), 
				row.names = FALSE, digits = 2, right = TRUE)
		
		cat("\n")
		#TODO classes
		#cat("\nClasses:\n", paste(rownames(x$beta[[1]]), sep = ", "))
		
	} else if(x$type == "cv") {
		
		cat("\nModels:\n\n")
		sel <- 1:5*floor(nmod(x)/5)
		
		err <- Err(x)
		feat <- colMeans(x$features)
		para <- colMeans(x$parameters)
		
		print(data.frame('Index: ' = sel, 
						'Lambda: ' = x$lambda[sel], 
						'Features: ' = print_with_metric_prefix(feat[sel]), 
						'Parameters: ' = print_with_metric_prefix(para[sel]), 
						'Error: ' = err[sel], check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)
		
		cat("\nBest model:\n\n")
		
		sel <- which(err == min(err))[1]
		
		print(data.frame('Index: ' = sel, 
						'Lambda: ' = x$lambda[sel], 
						'Features: ' = print_with_metric_prefix(feat[sel]), 
						'Parameters: ' = print_with_metric_prefix(para[sel]), 
						'Error: ' = err[sel], check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)
		
		cat("\n")
		
		#TODO classes
		#cat("\nClasses:\n", paste(unique(classes.true), sep = ", "))
		
	} else if(x$type == "subsampling") {
		
		cat("\nBest models:\n\n")
		
		err <- Err(x)
		
		if(nrow(err) <= 5) {
			sel <- 1:nrow(err)
		} else {
			sel <- 1:5*floor(ncol(err)/5)
		}
		
		model.sel <- apply(err, 1, function(y) which(min(y) == y)[1])
		feat <- sapply(1:length(sel), function(i) x$features[sel[i], model.sel[i]])
		para <- sapply(1:length(sel), function(i) x$parameters[sel[i], model.sel[i]])
		
		print(data.frame('Subsample: ' = sel, 
						'Model index: ' = model.sel,
						'Lambda: ' = x$lambda[model.sel], 
						'Features: ' = print_with_metric_prefix(feat), 
						'Parameters: ' = print_with_metric_prefix(para), 
						'Error: ' = sapply(1:length(sel), function(i) err[sel[i], model.sel[i]]), check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)
		
		cat("\n")
		
	} else if(x$type == "predict") {
		#TODO print some information
	
	} else {
		stop("Unknown type of sgl object")
	}
	
	
	#TODO msgl version check that  object$msgl_version exists
}


#' Print a numeric with metric prefix
#' 
#' @param x numeric to be printed 
#' @param digits number of significant digits 
#' 
#' @return a string 
#' 
#' @author Martin Vincent
#' @export
print_with_metric_prefix <- function(x, digits = 3) {
	
	if(length(x) > 1) {
		return(sapply(x, function(y) print_with_metric_prefix(y, digits = digits)))
	}
	
	metric_factor <-  c(1e+00, 1e+03, 1e+06, 1e+09)
	
	prefix <- c("", "k", "M", "G")

	sel <- max(which(x >= metric_factor))
	
	txt <- paste(round(x/metric_factor[sel], digits = digits), prefix[sel], sep="")
	
	return(txt)
}
