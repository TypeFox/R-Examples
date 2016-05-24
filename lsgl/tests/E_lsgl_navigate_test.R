#
#     Description of this R script:
#     R tests for linear multiple output sparse group lasso routines.
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


set.seed(100) #  ensures consistency of tests

## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 

N <- 50 #number of samples
p <- 50 #number of features
K <- 25  #number of groups

B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 

X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	

lambda<-lsgl.lambda(X,Y, alpha=1, lambda.min=.5, intercept=FALSE)

fit <-lsgl(X,Y, alpha=1, lambda = lambda, intercept=FALSE)

# print info
fit

# Test features
features(fit)

# parameters
parameters(fit)

nmod(fit)

models(fit)

coef(fit)