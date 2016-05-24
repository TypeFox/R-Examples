#
#     Description of this R script:
#     R test for linear multiple output using sparse group lasso routines.
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

library(lsgl)

# warnings = errors
options(warn=2)

set.seed(100) # This may be removed, it ensures consistency of the daily tests

## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 

N <- 50 #number of samples
p <- 25 #number of features
K <- 10  #number of groups

B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 
X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)

##Do cross validation
lambda <- lsgl.lambda(X1, Y1, alpha = 1, d = 25, lambda.min = 0.5, intercept = FALSE)

if(sgl.c.config()$omp.supported) {
	threads = 2L
} else {
	threads = 1L
}

fit.cv <- lsgl.cv(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE, max.threads = threads)

## Cross validation errors (estimated expected generalization error)
if(min(Err(fit.cv, loss = "SOVE")) > 0.05) stop()
