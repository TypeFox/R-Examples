#############################################################################
### Copyright Statement

### This R script is created by Meng Li 
### Email: mli9@ncsu.edu
#############################################################################


#############################################################################
### This file contains functions to generate skewed functional data with given 
### covariates and timepoints; mean, standard deviation and skewness functions;
### and correlation matrix of the latent Gaussian copula to model dependence. 
### The mean and standard deviation functions could be bivariate: depend both on 
### covariates and timepoints. 
#############################################################################




##   ORTHOGONAL BASIS SYSTEMS 
##   Legendre polynomials
## INPUT:  
##    t = the set of values, 
##		degree = the polynomial degree of the function
## 		normalized = logical value. 
##        If TRUE (default) then the values are normalized such that Phi^T Phi=I
## OUTPUT: 
##    vector of length $t$ 
##

legendre.polynomials=function(t, degree, normalized=TRUE){
  length.t= length(t)
  if (normalized==TRUE){
    if(degree==0){temp.lp=rep(1, length.t)}
    if(degree==1){temp.lp=sqrt(3)*(2*t-1)}
    if(degree==2){temp.lp=sqrt(5)*(6*t^2  - 6*t +1) }
    if(degree==3){temp.lp=sqrt(7)*(20*t^3 - 30*t^2  +12 *t -1) }
    temp.lp = temp.lp/sqrt(sum(temp.lp^2))
  }
  if (normalized!=TRUE){
    if(degree==0){temp.lp = rep(1, length.t)}
    if(degree==1){temp.lp =sqrt(3)*(2*t-1)}
    if(degree==2){temp.lp =sqrt(5)*(6*t^2  - 6*t +1) }
    if(degree==3){temp.lp =sqrt(7)*(20*t^3 - 30*t^2  +12 *t -1) }
  }
  return(temp.lp)
}


## 	 ORTHOGONAL BASIS SYSTEMS 
##   Discrete Fourier Transformation (DFT) basis system 
## INPUT:
##    t = the set of values,                       
##	  degree = the order of the function in the basis
## 		normalized = logical value. 
##      If  TRUE (=default) then the values are normalized such that Phi^T Phi=I
## OUTPUT: 
##    vector of length $t$ 
##
DFT.basis=function(t, degree=0, normalized=TRUE){
  length.t= length(t)
  parity.degree=degree - 2*floor(degree/2)
  degree.half = floor(degree/2)
  if (normalized==TRUE){
    if(degree==0){temp.lp=rep(1, length.t)}
    if((degree > 0)&&(parity.degree==0)){temp.lp=sqrt(2) * cos(2*degree.half*pi*t)}
    if((degree > 0)&&(parity.degree==1)){temp.lp=sqrt(2) * sin(2*(degree.half+1)*pi*t)}
    temp.lp = temp.lp/sqrt(sum(temp.lp^2))
  }
  if (normalized!=TRUE){
    if(degree==0){temp.lp=rep(1, length.t)}
    if((degree > 0)&&(parity.degree==0)){temp.lp=sqrt(2) * cos(2*degree.half*pi*t)}
    if((degree > 0)&&(parity.degree==1)){temp.lp=sqrt(2) * sin(2*(degree.half+1)*pi*t)}
  }
  return(temp.lp)
}

## Function: data.generator.y.F  
##    generates the process Y for the entire data set. 
##    require: DFT.basis , legendre.polynomials
##
## INPUT: 
##    n.subject = number of subjects in the data; 
##    n.timepoints = number of equally spaced locations in [0,1]
##	  s = mean vector of matrix;  
##    D = vector or matrix of log varaince 
##        D is a vector means the variance depends only on time; 
##        D is a matrix means the variance is bivariate wrt. time and covariate;
##    csi = vector or matrix for the shape; 
##		lambdas = vector of main eigenvalues; 
##    basis.system  = basys system used (default = Legendre)
##    var.noise = variance of white noise in the latent Gaussian Process (GP)
## OUTPUT: 
##    list(data = a n.subject x n.timepoints, matrix with one subject in each row; 
##         corr.true = the true corrlation matrix for the latent GP)
##

data.generator.y.F <- function(n.subject, n.timepoints, s , D, csi,	
                               lambdas = c(1/2, 1/4, 1/8), basis.system=legendre.polynomials, 
                               var.noise=0.10){
  
  # require(sn) #require the package sn
  
  # Generate the timepoints - equal-spaced in [0,1]
  temp.timepoints <- seq(0, 1, length = n.timepoints) 
  
  # Generate latent Gaussian Process
  sqrt.Lambda = diag(sqrt(lambdas))
  K.lambdas = length(lambdas)
  phi.basys.system = sapply(c(1:K.lambdas), function(k) basis.system(temp.timepoints, degree = k, normalized=FALSE))
  cov.true =  phi.basys.system %*% diag(lambdas)  %*% t(phi.basys.system)  +diag(var.noise, n.timepoints)      #noise here
  corr.true = cov2cor(cov.true)  # the true corrleation for latent process
  scores.matrix = matrix(rnorm(n.subject*K.lambdas), ncol=n.subject)
  var.pointwise = diag(phi.basys.system %*% diag(lambdas)  %*% t(phi.basys.system) )
  sigma = sqrt(var.pointwise+ var.noise)
  noise = sqrt(var.noise) * matrix(rnorm(n.subject*n.timepoints), nrow=n.subject)
  # lateng Gaussian Process : Z(t) in the paper
  Xd.data = t(phi.basys.system %*% sqrt.Lambda %*% scores.matrix ) + noise    
  if (n.subject==1)
    Xd.data = matrix(Xd.data , ncol=n.timepoints)       	
  
  # corpula: F_d(X_d)
  cdf.xi.data =  sapply(c(1:n.timepoints), function(k) pnorm(Xd.data[,k], mean=0, sd= sigma[k]))
  if (n.subject==1)
    cdf.xi.data = matrix(cdf.xi.data , ncol=n.timepoints)  	  
  
  # data set with location=0, scale=1, shape=csi 
  csi.ds <- csi
  if (length(csi.ds) == n.timepoints) {
    # when csi is a vector 
    yi.data <- sapply(c(1:n.timepoints), function(k) qsn(p=cdf.xi.data[,k], location=0, scale=1, shape=csi.ds[k], tol=1e-8)) 
    mean.yi.data <- as.vector(unlist(sapply(c(1:n.timepoints), function(k) sqrt(2/pi) * csi.ds[k] / sqrt( 1 + csi.ds[k]^2))))
    variance.yi.data <- 1 - mean.yi.data^2
    mean.matrix <- t(matrix(rep(mean.yi.data, n.subject), ncol=n.subject))
    sigma.matrix <-  t(matrix(rep(sqrt(variance.yi.data), n.subject), ncol=n.subject))
  } else {
    # when csi is a matrix
    yi.data <- sapply(1:length(cdf.xi.data), function(k) qsn(cdf.xi.data[k], location = 0, scale = 1, shape = csi.ds[k], tol = 1e-8))
    mean.matrix <- sqrt(2/pi) * csi.ds / sqrt( 1 + csi.ds^2 )
    variance.matrix <- 1 - mean.matrix^2
    sigma.matrix <- sqrt(variance.matrix)
  }

  if (length(s) == n.timepoints){ s = t(matrix(s, nrow=n.timepoints, ncol=n.subject))}
  if (length(D) == n.timepoints){ D = t(matrix(D, nrow=n.timepoints, ncol=n.subject))}

  # re-scale yi.data such that mean = s, logvariance = D
  final.result  =  s + (yi.data - mean.matrix)/sigma.matrix * exp(D/2)        

  list(data = final.result, corr.true = corr.true)
}

