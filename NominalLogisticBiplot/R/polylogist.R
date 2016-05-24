# file NominalLogisticBiplot/R/polylogist.R
# copyright (C) 2012-2013 J.L. Vicente-Villardon and J.C. Hernandez
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#Function that calculates an object with the fitted multinomial logistic regression for a nominal variable.
#----------------------Parameters--------------
  #y: response nominal variable.
  #x: matrix with independent variables.
  #penalization: value to correct the separation problem through the ridge regression
  #cte : it will be true if the model has an independent term.
  #tol : value to decide if the algorith should continue
  #maxiter : value to decide if the algorith should continue
  #show: boolean parameter if we want to see values in each iteration of the process.
polylogist <- function(y, x, penalization = 0.2, cte = TRUE, tol = 1e-04, maxiter = 200,show=FALSE) {

	# POLYLOGIST    Fits a polytomous ridge logistic regression
	# Produces a vector B of estimates of parameters
    # in a polytomoius logistic regression
    # The response vector must be constructed with integers starting at 0
    # the category labeled 0 is used for comparisons
    # The rows of B are for the categories and the column for the variables
    # To break the singularity of the information matrix a constant (penalization) is added
    # to the diagonal
  
  if(is.matrix(x)){
    n <- nrow(x)
  }else{
    n <- length(x)
  }	

	if (min(y)==1) y= y-1

	if (cte){                   
     x <- cbind(matrix(1,n),x)
  }
  
  p <- ncol(x)	

	J = max(y)
	Y = matrix(0, n, J)
	for (i in 1:n){ 
    if (y[i] > 0){
    		Y[i, y[i]] = 1
    }	   
	}

  beta <- matrix(0,J,p)

	yy = as.vector(Y)
	err = 0.1
	iter = 0
	m = matrix(0, n * J, n * J)
	prob = matrix(0, n, J)
	p0 = matrix(0, n, 1)
	
  xx = kronecker(diag(1, J), x)	   

	while ((err > tol) & (iter < maxiter)) {
		iter = iter + 1
		betaold = beta

		b = as.vector(t(beta))

		for (i in 1:n) {
			suma = 1
			for (j in 1:J) suma = suma + exp(sum(beta[j, ] * x[i, ]))
			for (j in 1:J) prob[i, j] = exp(sum(beta[j, ] * x[i, ]))/suma
			p0[i, 1] = 1/suma
		}

		P = as.vector(prob)

		for (j in 1:J)
     for (k in 1:J){
    			mj = matrix(0, n, n)
    			for (i in 1:n) {
    				if (k == j) {
    					mj[i, i] = prob[i, j] * (1 - prob[i, j])
    				}
    				if (!(k == j)) {
    					mj[i, i] = -1 * prob[i, j] * prob[i, k]
    				}
    			}
    			m[((j - 1) * n + 1):(j * n), ((k - 1) * n + 1):(k * n)] = mj
		 }
		a = (t(xx) %*% m %*% xx) + 2 * diag(J*p) * penalization
		u = (t(xx) %*% matrix((yy - P), n * J, 1)) - 2 * b * penalization
		b = b + solve(a) %*% u
		beta = t(matrix(b, p, J))
		err = sum(sum((betaold - beta) * (betaold - beta)))
		#if (show) print(c(iter, err))
	}


	y0=matrix(0,n,1)
	for (i in 1:n) if (sum(Y[i,]) == 0) y0[i]=1
  cov=solve(a)
  model<-list()
  model$fitted=cbind(p0,prob)
  model$cov=cov

  model$Y=cbind(y0,Y)
  model$beta=beta
  model$stderr=sqrt(t(matrix(diag(cov), p, J)))
  model$logLik=sum(log(model$fitted^model$Y))
  model$Deviance = -2*sum(log(model$fitted^model$Y))
  model$AIC = model$Deviance + 2*nrow(beta)*ncol(beta)
  model$BIC = model$Deviance + nrow(beta)*ncol(beta)*log(n)

  class(model)='polylogist'
  return(model)
}
