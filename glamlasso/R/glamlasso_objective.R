#
#     Description of this R script:
#     R interface for the function objective in the glamlasso package.
#
#     Intended for use with R.
#     Copyright (C) 2015 Adam Lund
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
 
#' @name objective
#' 
#' @aliases glamlasso_objective 
#' 
#' @title Compute objective values 
#' 
#' @description Computes the objective values of the lasso penalized log-likelihood problem
#'              for the models implemented in the package glamlasso.
#'              
#' @usage objective(Y, 
#'           Weights, 
#'           X1, X2, X3, 
#'           Beta, 
#'           lambda,
#'           penalty.factor, 
#'           family)
#' 
#' @param Y The response values, a 3d array of size \eqn{n_1 \times n_2 \times n_3}. 
#' @param Weights Observation weights, a 3d array of size \eqn{n_1 \times n_2 \times n_3}.
#' @param X1,X2,X3 The three tensor components of the tensor design matrix, each of size
#'  \eqn{n_i \times p_i}, \eqn{i = 1, 2, 3}.
#' @param Beta A coefficient matrix of size \eqn{p_1p_2p_3 \times }\code{nlambda}.
#' @param lambda The sequence of penalty parameters for the regularization path.
#' @param penalty.factor A 3d array of size  \eqn{p_1 \times p_2 \times p_3}. Is multiplied with each 
#' element in  \code{lambda} to allow differential shrinkage on the coefficients. 
#Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. 
#Default is 1 for all variables (and implicitly infinity for variables listed in exclude). Note: the penalty 
#factors are internally rescaled to sum to nvars, and the lambda sequence will reflect this change.
#' @param family A string indicating the model family (essentially the response distribution).
#' 
#' @return
#' A vector of length \code{length(lambda)} containing the objective values for each \code{lambda} value. 
#' 
#' @examples
#' \dontrun{
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 13; p2 <- 5; p3 <- 4
#' X1 <- matrix(rnorm(n1 * p1), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3), n3, p3) 
#' Beta <- array(rnorm(p1 * p2 * p3) * rbinom(p1 * p2 * p3, 1, 0.1), c(p1 , p2, p3))
#' mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rnorm(n1 * n2 * n3, mu), dim = c(n1, n2, n3))
#' fit <- glamlasso(X1, X2, X3, Y, family = "gaussian", iwls = "exact")
#' objfit <- objective(Y, NULL, X1, X2, X3, fit$coef, fit$lambda, NULL, fit$family)
#' plot(objfit, type = "l")
#' }

objective <-function(Y, Weights, X1, X2, X3, Beta, lambda, penalty.factor, family) {
  
##get dimensions of problem
dimX <- rbind(dim(X1), dim(X2), dim(X3))
n <- prod(dimX[ , 1])
p <- prod(dimX[ , 2])
n1 <- dimX[1, 1]
n2 <- dimX[2, 1]
n3 <- dimX[3, 1]
p1 <- dimX[1, 2]
p2 <- dimX[2, 2]
p3 <- dimX[3, 2]

##check if lambda is specified
if(is.null(lambda)){stop(paste("no lambda sequence is specified"))}

##check if penalty.factor is specified; if not set to all ones
if(is.null(penalty.factor)){
  
penalty.factor <- matrix(1, p1, p2 * p3)

}else if(prod(dim(penalty.factor)) != p){
  
stop(
paste("number of elements in penalty.factor (", length(penalty.factor),") is not equal to the number of coefficients (", p,")", sep = "")
)

}else {penalty.factor <- matrix(penalty.factor, p1, p2 * p3)}  

##check if weights/nof trials is specfied for the binomial model
if(family == "binomial" & is.null(Weights)){

stop(paste("for binomial model number of trials (Weights) must be specified"))

}

if(is.null(Weights)){
  
Weights <- matrix(1, n1, n2 * n3)

}else{Weights <- matrix(Weights, n1, n2 * n3)}

Y <- matrix(Y, n1, n2 * n3)
BetaArr <- array(Beta, c(p1, p2 * p3, length(lambda)))

##get objective values
out <- getobj(Y, Weights,  X1, X2, X3, BetaArr, lambda, penalty.factor, family)$Obj

return(out)

}