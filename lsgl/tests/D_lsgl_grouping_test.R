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

N <- 25 #number of samples
p <- 25 #number of features
K <- 10  #number of groups

B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K) 

X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)	

grouping <- rep(LETTERS[1:5],5)

lambda<-lsgl.lambda(X,Y, grouping = grouping, alpha=0, lambda.min=0.1, intercept=FALSE)

fit <-lsgl(X,Y, grouping = grouping, alpha=0, lambda = lambda, intercept=FALSE)

if(min(Err(fit, X)) > 1) stop()

tmp <- which(rowSums(abs(fit$beta[[2]])) != 0)
if(! all((tmp[1]+5*1:4) %in% tmp)) stop()

## Test single fit i.e. K = 1
y <- Y[,1]

lambda<-lsgl.lambda(X,y, grouping = grouping,alpha=0, lambda.min=.5, intercept=FALSE)
fit <-lsgl(X, y, grouping = grouping, alpha=0, lambda = lambda, intercept=FALSE)
res <- predict(fit, X)

tmp <- which(rowSums(abs(fit$beta[[2]])) != 0)
if(! all((tmp[1]+5*1:4) %in% tmp)) stop()

