#
#     Description of this R script:
#     R interface for linear multi-response sparse group lasso routines.
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

#' @title Predict
#'
#' @description 
#' Compute the predicted response matrix for a new data set.
#'
#' @param object an object of class lsgl, produced with \code{lsgl}.
#' @param x a data matrix of size \eqn{N_\textrm{new} \times p}.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param ... ignored.
#' @return
#' \item{Yhat}{the predicted response matrix (of size \eqn{N_\textrm{new} \times K})}
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
#' X2<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y2 <-X2%*%B+matrix(rnorm(N*K,0,1),N,K)
#' 
#' #### Fit models using X1
#' lambda <- lsgl.lambda(X1, Y1, alpha = 1, d = 25L, lambda.min = 0.5, intercept = FALSE)
#' fit <- lsgl(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#' 
#' # Predict Y2 using the estimated models and X2
#' res <- predict(fit, X2)
#' 
#'  # Frobenius norm \|Y2hat - Y2\|_F
#' sapply(res$Yhat, function(y) sqrt(sum((y - Y2)^2)))
#'
#' # This is proportional to the errors compute with Err:
#' Err(fit, X2, Y2)*length(Y2)
#'
#' @author Martin Vincent
#' @method predict lsgl
#' @importFrom methods is
#' @importFrom methods as
#' @export
#' @useDynLib lsgl, .registration=TRUE
predict.lsgl <- function(object, x, sparse.data = is(x, "sparseMatrix"), ...) 
{
	# Get call
	cl <- match.call()
	
	if(is.null(object$beta)) stop("No models found -- missing beta")
	
	if(object$intercept){
		# add intercept
		x <- cBind(Intercept = rep(1, nrow(x)), x)
	}	
	
	#Check dimension of x
	if(dim(object$beta[[1]])[1] != ncol(x)) stop("x has wrong dimension")
	
	object$beta <- lapply(object$beta, t)
	
	data <- list()
	data$sparseX <- FALSE
	
	if(is(x, "kron")) {
		
		if(length(x) == 2) {
			callsym <- "lsgl_kdx"
		} else if(length(x) == 3) {
			callsym <- "lsgl_ktx"
		} else {
			stop("unsupported kron")
		}
		
		data$X <- x
		res <- sgl_predict(callsym, "lsgl", object, data)
				
	} else if(sparse.data) {
		
		data$X <- as(x, "CsparseMatrix")
		data$sparseX <- TRUE
		
		res <- sgl_predict("lsgl_xs_yd", "lsgl", object, data)
		
	} else {
		
		data$X <- as.matrix(x)
		
		res <- sgl_predict("lsgl_xd_yd", "lsgl", object, data)
		
	}
	
	#Responses
	
	res$Yhat <- lapply(res$responses$link, t)
	res$responses <- NULL
	res$call <- cl
	
	return(res)
}
