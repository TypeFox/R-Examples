##
##  PURPOSE:   Truncated multivariate normal distribution
##             * random numbers generation via Gibbs sampling
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   20/11/2007
##
##  FUNCTION:  rTMVN
##
## ======================================================================


## *************************************************************
## rTMVN
## *************************************************************
rTMVN <- function(n, mean=c(0, 0), Sigma=diag(2), a, b, trunc, xinit)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")

  nx <- length(mean)
  if (nx == 1){    
    if (length(Sigma) != 1) stop(paste("length(mean)=", nx, ", but length(Sigma)=", length(Sigma), "!", sep=""))
    if (missing(trunc)){
      trunc <- 4
      a <- b <- 0
    }else{
      if (trunc == 0) b <- 0
      else if (trunc == 1) b <- 0
           else if (trunc == 2) b <- 0
                else if (trunc == 4) a <- b <- 0      
      if (missing(a) | missing(b)) stop("a and b must be given")      
    }
    
    return(rTNorm(n=n, mean=mean, sd=sqrt(Sigma), a=a, b=b, trunc=trunc))    
  }
  
  if (is.null(dim(Sigma))) stop(paste("length(mean)=", nx, ", but dim(Sigma) is NULL!", sep=""))
  Q <- chol2inv(chol(Sigma))       ## To check whether Sigma is PD

  if (missing(trunc)){
    trunc <- rep(4, nx)
    a <- b <- rep(0, nx)
  }else{
    if (length(trunc) != nx) stop(paste("length(mean)=", nx, ", but length(trunc)=", length(trunc), "!", sep=""))
    if (any(!(trunc %in% c(0, 1, 2, 3, 4)))) stop("All trunc values must be from {0, 1, 2, 3, 4}.")        
    if (!any(trunc==3)) b <- rep(0, nx)
    if (sum(trunc==4) == nx) a <- b <- rep(0, nx)        
    if (missing(a) | missing(b)) stop("a and b must be given")
    if (length(trunc) != length(a) | length(trunc) != length(b)) stop("a, b and trunc must have the same length")
    a[trunc==4] <- 0
    b[trunc==0] <- 0
    b[trunc==1] <- 0    
    b[trunc==2] <- 0
    b[trunc==4] <- 0    
    if (any(a[trunc==3] >= b[trunc==3])) stop("a must be lower than b when trunc = 3")     
  }    

  if (missing(xinit)){
    sigmas <- sqrt(diag(Sigma))
    xinit <- mean
    xinit[trunc==0] <- a[trunc==0] + sigmas[trunc==0]
    xinit[trunc==1] <- a[trunc==1]
    xinit[trunc==2] <- a[trunc==2] - sigmas[trunc==2]    
    xinit[trunc==3] <- (a[trunc==3] + b[trunc==3])/2
  }else{
    if (length(xinit) != nx) stop(paste("length(mean)=", nx, ", but length(xinit)=", length(xinit), "!", sep=""))
    if (any(xinit[trunc==0] <= a[trunc==0])) stop("xinit must be higher than a when trunc = 0")
    if (any(xinit[trunc==1] != a[trunc==1])) stop("xinit must be equal to a when trunc = 1")    
    if (any(xinit[trunc==2] >= a[trunc==2])) stop("xinit must be lower than a when trunc = 2")                     
    if (any(xinit[trunc==3] >= b[trunc==3])) stop("xinit must be lower than b when trunc = 3")
    if (any(xinit[trunc==3] <= a[trunc==3])) stop("xinit must be higher than a when trunc = 3")             
  }  

  Sigma <- Sigma[lower.tri(Sigma, diag=TRUE)]

  SAMPLE <- .C("rTMVN1_R", x      =double(nx*n),
                           beta   =double(nx),
                           sigmaR2=double(nx),
                           L      =double((nx-1)*nx/2),
                           err    =integer(1),
                           xinit  =as.double(xinit),
                           mu     =as.double(mean),
                           Sigma  =as.double(Sigma),
                           a      =as.double(a),
                           b      =as.double(b),
                           trunc  =as.integer(trunc),
                           p      =as.integer(nx),
                           npoints=as.integer(n),
               PACKAGE=thispackage)

  SAMPLE$x <- matrix(SAMPLE$x, ncol=nx, nrow=n, byrow=TRUE)
  colnames(SAMPLE$x) <- names(mean)  
  
  return(SAMPLE$x)
}  


