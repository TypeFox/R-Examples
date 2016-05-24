######################################################
### gaussDiv.R
### implement multiple divergence measures to compare normal
### pdfs, i.e. similarity and dissimilarity measures
### Henning Rust                        Paris 18/02/09

### all the implemented similarity and dissimilarity measures
### are described in the chapter
### Dissimilatiry Measures for Probability Distributions
### in the book "Analysis of Symbolic Data" by Hans-Hermann Bock 

## visible functions
##--------------------------

### general wrapper function
normdiff <- function(mu1,sigma1=NULL,mu2,sigma2=sigma1,inv=FALSE,s=0.5,method=c("Mahalanobis","KL","J","Chisq","Hellinger","L2","Euclidean")){
  ## maybe I catch some erros, checking simgas for square and symmetry
  ## Euclidean for having this also in the same wrapper
  d <- switch(match.arg(method),
              Mahalanobis=normdiff.maha(mu1,sigma1,mu2),
              KL=normdiff.KL(mu1,sigma1,mu2,sigma2),
              J=normdiff.J(mu1,sigma1,mu2,sigma2),
              Chisq=normdiff.Chisq(mu1,sigma1,mu2,sigma2),
              Hellinger=normdiff.Hellinger(mu1,sigma1,mu2,sigma2),
              L2=normdiff.L2(mu1,sigma1,mu2,sigma2),
              Euclidean=sqrt(sum((mu1-mu2)**2)))
  class(d) <- "normdiff"
  attr(d,"method") <- match.arg(method)
  if(!is.null(inv)) attr(d,"inverse") <- inv
  if(!is.null(s)) attr(d,"s") <- s
  return(d)
}


## internal functions
##--------------------------

### trace of a matrix
tt <- function(A) sum(diag(A))

### Mahalanobis distance
maha <- function(x,A,factor=FALSE){
  ## if factor, then distance is multiplied with 0.5
  inv.A <- solve(A)
  d <- t(x)%*%inv.A%*%x
  if(factor) d <- 0.5*d
  return(as.vector(d))
}

### Mahalanobis Distance
normdiff.maha <- function(mu1,sigma1,mu2){
  d <- NA
  if(!(is.null(sigma1)))
    d <- maha(mu1-mu2,sigma1,factor=TRUE)
  return(as.vector(d))
}

### Kullback-Leibler divergence
normdiff.KL <- function(mu1,sigma1,mu2,sigma2){
  d <- NA
  if(!(is.null(sigma1)|is.null(sigma2))){
    N <- nrow(sigma1)
    inv.sigma1 <- solve(sigma1)
    d <- maha(mu1-mu2,sigma1)+
      tt(inv.sigma1%*%sigma2-diag(1,N)) +
        log(det(sigma1)/det(sigma2))
    d <- 0.5*d
  }
  return(as.vector(d))
}

### J-coefficient, symmetric
### symmetrized Kullback-Leibler
normdiff.J <- function(mu1,sigma1,mu2,sigma2){
  d <- NA
  if(!(is.null(sigma1)|is.null(sigma2))){
    N <- nrow(sigma1)
    inv.sigma1 <- solve(sigma1)
    inv.sigma2 <- solve(sigma2)
    d <- maha(mu1-mu2,sigma1)+maha(mu1-mu2,sigma2)+
      tt(inv.sigma1%*%sigma2-diag(1,N)) +
        tt(inv.sigma2%*%sigma1-diag(1,N))
    d <- 0.5*d
  }
  return(as.vector(d))
}

### chi-square divergence for pdfs
normdiff.Chisq <- function(mu1,sigma1,mu2,sigma2){
  d <- NA
  if(!(is.null(sigma1)|is.null(sigma2))){
    N <- nrow(sigma1)
    inv.sigma1 <- solve(sigma1)
    inv.sigma2 <- solve(sigma2)
    sig1.invsig2 <- sigma1%*%inv.sigma2
    d <- det(sig1.invsig2)/sqrt(det(2*sig1.invsig2-diag(1,N)))*
      exp(0.5*(maha(2*inv.sigma2%*%mu2-inv.sigma1%*%mu1,2*inv.sigma2-inv.sigma1)+
               maha(mu1,sigma1)-2*maha(mu2,sigma2)))-1
  }
  return(as.vector(d))
}

### Hellinger distance,
### similarity measure
### 0<d<1,
### d=1, if P=Q,
### d=0, if P, Q have disjoint supports
### symmetric for s=0.5
normdiff.Hellinger<- function(mu1,sigma1,mu2,sigma2,s=0.5,inv=FALSE){
  d <- NA
  if(!(is.null(sigma1)|is.null(sigma2))){
    N <- nrow(sigma1)
    I <- diag(1,N)
    inv.sigma1 <- solve(sigma1)
    inv.sigma2 <- solve(sigma2)
    sig1inv.sig2 <- inv.sigma1%*%sigma2
    sig2inv.sig1 <- inv.sigma2%*%sigma1
    d <- det(s*I+(1-s)*sig1inv.sig2)**(-s/2)*
      det((1-s)*I+s*sig2inv.sig1)**(-(1-s)/2)*
        exp(0.5*(maha(s*inv.sigma2%*%mu2+(1-s)*inv.sigma1%*%mu1,s*inv.sigma2+(1-s)*inv.sigma1)-
                 s*maha(mu2,sigma2)-(1-s)*maha(mu1,sigma1)))
    if(inv) d <- 1-d
  }
  return(as.vector(d))
}

### Minkowskis L2-distance
### symmetric dissimilarity coefficient
normdiff.L2<- function(mu1,sigma1,mu2,sigma2){
  d <- NA
  if(!(is.null(sigma1)|is.null(sigma2))){
    N <- nrow(sigma1)
    inv.sigma1 <- solve(sigma1)
    inv.sigma2 <- solve(sigma2)
    d <- 1/(2**N*pi**(N/2))*(1/sqrt(det(sigma1))+1/sqrt(det(sigma2)))-
      2/((2*pi)**(N/2)*sqrt(det(sigma1+sigma2)))*
        exp(0.5*(maha(inv.sigma1%*%mu1+inv.sigma2%*%mu2,inv.sigma1+inv.sigma2)-
                 maha(mu1,sigma1)-maha(mu2,sigma2)))
  }
  return(as.vector(d))
}


### a print function
print.normdiff <- function(x,...){
  cat("Gaussian PDF distance:",attr(x,"method"))
  if(attr(x,"method")=="Hellinger"){
    cat(" with s=",attr(x,"s"),", ")
    if(attr(x,"inv")==TRUE)
      cat("inverted")
  }
  cat("\n")
  cat(x,"\n")
}
