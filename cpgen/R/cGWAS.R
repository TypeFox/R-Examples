#
# cGWAS.R
# Claas Heuer, June 2014
#
# Copyright (C)  2014 Claas Heuer
#
# This file is part of cpgen.
#
# cpgen is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# cpgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# in file.path(R.home("share"), "licenses").  If not, see
# <http://www.gnu.org/licenses/>.
#

# cGWAS

cGWAS <- function(y,M,X=NULL,V=NULL,dom=FALSE, verbose=TRUE){

if(!is.vector(y) | !is.numeric(y)) stop("y must be a numeric vector")
y <- as.numeric(y)
id = 1:length(y)
isy = id[!is.na(y)]
n = length(isy)
if(n<length(y)) { stop("NA's in y are not allowed")}
sparse=FALSE

if(is.null(X)) X<-array(1,dim=c(n,1))
if(is.null(V)) V<-sparseMatrix(i=1:n,j=1:n,x=rep(1,n))
if(class(X)!="matrix") stop("X must be of type 'matrix'")
if(anyNA(X)) stop("NAs in X are not allowed")
if(class(M)!="matrix") stop("M must be of type 'matrix'")
if(anyNA(M)) stop("NAs in M are not allowed")

if (class(V) == "dgCMatrix"){
  sparse=TRUE
} else { 
    if (class(V) != "matrix") { stop("V must be either of type 'matrix' or 'dgCMatrix'") } 
  }


## this is only for internal use, namely: cGWAS.emmax
second_transform=FALSE;
V2 = array(1,dim=c(1,1))  	 

#	 .Call( "cGWAS", y[isy], M[isy,], V[isy,isy], X[isy,], dom, sparse, threads, PACKAGE = "cpgen" )
gwa <- .Call( "cGWAS", y, M, V, V2, X, dom, second_transform, sparse, verbose, options()$cpgen.threads, PACKAGE = "cpgen" )

if(dom) { 

  colnames(gwa$p_value) <- c("add","dom")
  colnames(gwa$beta) <- c("add","dom")
  colnames(gwa$se) <- c("add","dom") 

} else {

  gwa$p_value <- gwa$p_value[,1]
  gwa$beta <- gwa$beta[,1]
  gwa$se <- gwa$se[,1]

  }

return(gwa)

}



# cGWAS.emmax

cGWAS.emmax <- function(y,M,A=NULL,X=NULL,dom=FALSE,verbose=TRUE,scale_a = 0, df_a = -2, scale_e = 0, df_e = -2,niter=15000,burnin=7500,seed=NULL){

if(!is.vector(y) | !is.numeric(y)) stop("y must be a numeric vector")
y <- as.numeric(y)
id = 1:length(y)
isy = id[!is.na(y)]
n = length(isy)
if(n<length(y))   stop("NAs in y are not allowed")
if(is.null(X)) X <- array(1,dim=c(n,1))

if(class(X)!="matrix") stop("X must be of type 'matrix'")
if(anyNA(X)) stop("NAs in X are not allowed")
if(class(M)!="matrix") stop("M must be of type 'matrix'")
if(anyNA(M)) stop("NAs in M are not allowed")

if(verbose) cat("\nComputing Eigen Decomposition and Estimating V\n")

lambda=0.01
if(is.null(A)) A = cgrm(M,lambda=lambda)
if(class(A)!="matrix") stop("A must be of type 'matrix'")
if(anyNA(A)) stop("NAs in A are not allowed")

UD <- eigen(A)

UX <- t(UD$vectors) %c% X
Uy <- (t(UD$vectors) %c% y)[,1]
D_sqrt <- sqrt(UD$values)
Z <- sparseMatrix(i=1:n,j=1:n,x=D_sqrt)

par_random <- list(list(scale=scale_a,df=df_a,sparse_or_dense="sparse",method="ridge"))
if(missing(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }

# set the number of threads to 1 for clmm
old_threads <- get_num_threads()
set_num_threads(1,silent=TRUE)

mod <- clmm(y = Uy, X = UX, Z = list(Z), par_random = par_random, scale_e=scale_e, df_e=df_e, niter=niter, burnin=burnin, seed=seed, verbose=FALSE)

# set number of threads to old value for GWAS
set_num_threads(old_threads,silent=TRUE)

s2a = mod[[4]]$posterior$variance_mean
s2e = mod$Residual_Variance$Posterior_Mean
if(verbose) cat(paste("\nVariance Components:\nMarker: ",round(s2a,digits=3),"\nResidual: ",round(s2e,digits=3),"\n",sep=""))

v = 1/sqrt(UD$values*s2a + s2e)
V = sparseMatrix(i=1:n,j=1:n,x=v)
sparse=TRUE
if(verbose) cat("\nRunning GWAS\n")
second_transform=TRUE

gwa<- .Call( "cGWAS", Uy, M, V, t(UD$vectors), UX, dom, second_transform, sparse,verbose, options()$cpgen.threads, PACKAGE = "cpgen" )

if(dom) { 

  colnames(gwa$p_value) <- c("add","dom")
  colnames(gwa$beta) <- c("add","dom")
  colnames(gwa$se) <- c("add","dom") 

} else {

  gwa$p_value <- gwa$p_value[,1]
  gwa$beta <- gwa$beta[,1]
  gwa$se <- gwa$se[,1]

  }

gwa$marker_variance = s2a
gwa$residual_variance = s2e

return(gwa)

}




