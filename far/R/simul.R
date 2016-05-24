#  *****************************************************************************
#   File : simul.R
#         ************************************************************
#   Description : 
#       Simulation of functional processes.
#       FAR, FARX, Wiener,...
#   Version : 1.3
#   Date : 2005-01-11
#         ************************************************************
#   Author : Julien Damon <julien.damon@gmail.com>
#   License : LGPL
#   URL: https://github.com/Looping027/far
#  *****************************************************************************

#  *****************************************************************************
#   Title : interpol.matrix
#         ************************************************************
#   Description : 
#       Calculate the matrix giving the linear interpolation of values
#       of a function
#   Version : 0.3
#   Date : 2003-04-06
#  *****************************************************************************
interpol.matrix <- function(n=12,m=24,tol=sqrt(.Machine$double.eps))
{
  # Inner function to calculate the interpolation
  .inner.interpol.matrix <- function(n,m,tol=sqrt(.Machine$double.eps))
    {
      if ((n<=0) || (m<=0)) return(NULL)
      temp <- outer((1+0:(n-1))/n,(1+0:(m-1))/m,"-")
      temp <- m*(1/m-abs(temp))
      temp[temp<tol] <- 0
      if (n>m) temp[,m] <- temp[,m]/apply(temp,1,"sum")
      return(temp)
    }

  # Handling errors
  n[n<0] <- 0
  m[m<0] <- 0
  size <- max(length(n),length(m))
  n <- c(n,rep(0,size-length(n)))
  m <- c(m,rep(0,size-length(m)))

  # Main calculation
  if (size==0) return(NULL)
  if (size==1) {
    result <- .inner.interpol.matrix(n,m,tol)
  } else {
    n2 <- cumsum(n)
    m2 <- cumsum(m)
    result <- matrix(0,nrow=n2[size],ncol=m2[size])
    n2 <- c(0,n2)
    m2 <- c(0,m2)
    for (k in 1:size){
      temp <- .inner.interpol.matrix(n[k],m[k],tol)
      if (!is.null(temp))
      result[(n2[k]+1):n2[k+1],(m2[k]+1):m2[k+1]] <- temp
    }
  }

  # returning result
  return(result)
}

#  *****************************************************************************
#   Title : BaseK2BaseC
#         ************************************************************
#   Description : 
#       Given the coordinates in the Karhunen Loeve expansion
#       base of the Wiener, compute the coordinates in the
#       canonical base
#   Version : 1.0
#   Date : 2001-03-16
#  *****************************************************************************
BaseK2BaseC <- function(x,nb=nrow(x))
{
    if (!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    cst1 <- ((1:n)-0.5)*pi
    prod1 <- outer((1:nb)/nb,cst1,"*")
    res <- sin(prod1) %*% (x / cst1) * sqrt(2)
    res <- as.fdata(res,dates=1:ncol(x))
    return(res)
}

#  *****************************************************************************
#   Title : simul.wiener
#         ************************************************************
#   Description : 
#       Wiener simulation using the Karhunen Loeve expansion
#   Version : 1.1
#   Date : 2003-04-06
#  *****************************************************************************
simul.wiener <- function(m=64,n=1,m2=NULL) {
    if (is.null(m2)) m2 <- 2*m
    noise <- matrix(rnorm(n*m2),nrow=m2,ncol=n)
    return(BaseK2BaseC(noise,nb=m))
}

#  *****************************************************************************
#   Title : simul.far.wiener
#         ************************************************************
#   Description : 
#       FAR(1) simulation using the Karhunen Loeve expansion
#       of the Wiener noise
#   Version : 1.2
#   Date : 2003-06-13
#  *****************************************************************************
simul.far.wiener<-function(m=64,n=128,
            d.rho=diag(c(0.45,0.90,0.34,0.45)),cst1=0.05,m2=NULL)
{
    if (is.null(m2)) m2 <- 2*m
    if (ncol(d.rho) > m2) d.rho <- d.rho[,1:m2,drop=FALSE]
    if (nrow(d.rho) > m2) d.rho <- d.rho[1:m2,,drop=FALSE]    
    if (ncol(d.rho) < m2)
    {
        d2.rho <- diag(cst1/((1:m2)^2)+(1-cst1)/exp(1:m2))
        d2.rho[1:nrow(d.rho),1:ncol(d.rho)] <- d.rho
        d.rho <- d2.rho
    }
    noise <- matrix(rnorm(n*m2),nrow=m2,ncol=n)
    x <- matrix(0,nrow=m2,ncol=n)
    x[,1] <- noise[,1]
    for (i in (2:n))
        x[,i] <- d.rho %*% x[,i-1] + noise[,i-1]
    return(BaseK2BaseC(x,nb=m))
}

#  *****************************************************************************
#   Title : simul.far.sde
#         ************************************************************
#   Description : 
#       Far(1) simulation using a stochastic diffential equation
#   Version : 1.0
#   Date : 2001-03-16
#  *****************************************************************************
simul.far.sde <- function(coef=c(0.4,0.8),n=80,p=32,sigma=1)
{
    m <- 10
    deriv <- rep(0,(n+m)*p)
    x <- rep(0,(n+m)*p)
    delta <- 1/p
    noise <- sigma * sqrt(delta) * rnorm((n+m)*p)
    for (i in 2:((n+m)*p))
    {
        deriv[i] <- deriv[i-1]  + noise[i] -
                (coef[2] * deriv[i-1] + coef[1] * x[i-1]) * delta
        x[i] <- x[i-1] + delta * deriv[i-1]
    }
    dim(x) <- c(p,n+m)
    x <- x[,(m+1):(n+m),drop=FALSE]
    x <- as.fdata(x,dates=1:ncol(x))
    return(x)
}

#  *****************************************************************************
#   Title : simul.far
#         ************************************************************
#   Description : 
#       FAR(1) simulation
#   Version : 1.4
#   Date : 2003-06-13
#  *****************************************************************************
simul.far <- function(m=12,n=100,base=base.simul.far(24,5),
             d.rho=diag(c(0.45,0.90,0.34,0.45)),
             alpha=diag(c(0.5,0.23,0.018)),
             cst1=0.05)
{
    m2 <- nrow(base)

    if (ncol(d.rho) > m2) d.rho <- d.rho[,1:m2,drop=FALSE]
    if (nrow(d.rho) > m2) d.rho <- d.rho[1:m2,,drop=FALSE]    

    if (ncol(alpha) > m2) alpha <- alpha[,1:m2,drop=FALSE]
    if (nrow(alpha) > m2) alpha <- alpha[1:m2,,drop=FALSE]    

    # Initialization
    # X associated projector computation
    ProjX <- interpol.matrix(m,m2) %*% orthonormalization(base)
    
    # Matrix of the linear form
    mat.rho <- diag(cst1/((1:m2)^2)+(1-cst1)/exp(1:m2))
    mat.rho[1:nrow(d.rho),1:ncol(d.rho)] <- d.rho
    
    # Coefficients
    Alpha <- diag(cst1/(1:m2))
    Alpha[1:nrow(alpha),1:ncol(alpha)] <- m2*alpha
    Tmpmat <- Alpha - (mat.rho %*% Alpha %*% t(mat.rho))
    Tmpmat <- Tmpmat[1:nrow(alpha),1:nrow(alpha),drop=FALSE]
    test <- (try(Tmpmat <- t(pdMatrix(pdSymm(Tmpmat),factor=TRUE)),TRUE))
    if (inherits(test,"try-error"))
      stop("The resulting matrix of the noise is not positive definite.\n Modify d.rho and/or alpha\n to get (alpha - (rho %*% alpha %*% t(rho))) positive definite.")
    Mu <- diag(cst1/(1:m2))
    Mu[1:nrow(alpha),1:nrow(alpha)] <- Tmpmat

    # Strong white noise simulation
    epsilon <- matrix(rnorm(m2*n),nrow=m2,ncol=n)
    epsilon <- Mu %*% epsilon
    
    # Creation of an empty X
    X <- matrix(0,nrow=m2,ncol=n)
    
    # Iteration of simulation
    X[,1] <- epsilon[,1]
    for (i in (2:n))
        X[,i] <- mat.rho %*% X[,i-1] + epsilon[,i]
    
    # Projection in the canonical basis
    as.fdata(ProjX %*% X,dates=1:n)
}

#  *****************************************************************************
#   Title : base.simul.far
#         ************************************************************
#   Description : 
#       Example of base used to simulate FAR(1)
#   Version : 1.0
#   Date : 2001-07-02
#  *****************************************************************************
base.simul.far <- function(m=24,n=5)
{
    res <- matrix(0,ncol=n,nrow=m)
    xval <- (0:(m-1)) / m
    for (i in 1:n)
        res[,i] <- sin(pi*i*xval)
    res
}

#  *****************************************************************************
#   Title : simul.farx
#         ************************************************************
#   Description : 
#       FARX(1,1) simulation using a strong white noise
#   Version : 1.4
#   Date : 2003-06-13
#  *****************************************************************************
simul.farx <- function(m=12,n=100,base=base.simul.far(24,5),
                base.exo=base.simul.far(24,5),
                d.a=matrix(c(0.5,0),nrow=1,ncol=2),
                alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
                d.rho=diag(c(0.45,0.90,0.34,0.45)),
                alpha=diag(c(0.5,0.23,0.018)),
                d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
                cst1=0.05)
{
    if (length(m)==1) m <- rep(m,2)
    m2 <- c(nrow(base),nrow(base.exo))

    if (ncol(d.a) > m2[2]) d.a <- d.a[,1:m2[2],drop=FALSE]
    if (nrow(d.a) > m2[1]) d.a <- d.a[1:m2[1],,drop=FALSE]    

    if (ncol(d.rho) > m2[1]) d.rho <- d.rho[,1:m2[1],drop=FALSE]
    if (nrow(d.rho) > m2[1]) d.rho <- d.rho[1:m2[1],,drop=FALSE]    

    if (ncol(alpha) > m2[1]) alpha <- alpha[,1:m2[1],drop=FALSE]
    if (nrow(alpha) > m2[1]) alpha <- alpha[1:m2[1],,drop=FALSE]    

    if (ncol(alpha.conj) > m2[2]) alpha.conj <- alpha.conj[,1:m2[2],drop=FALSE]
    if (nrow(alpha.conj) > m2[1]) alpha.conj <- alpha.conj[1:m2[1],,drop=FALSE]    

    if (ncol(d.rho.exo) > m2[2]) d.rho.exo <- d.rho.exo[,1:m2[2],drop=FALSE]
    if (nrow(d.rho.exo) > m2[2]) d.rho.exo <- d.rho.exo[1:m2[2],,drop=FALSE]    

    # Creation of the global basis
    basetot <- matrix(0,ncol=sum(m2),nrow=sum(m2))
    basetot[1:m2[1],1:m2[1]] <- orthonormalization(base)
    basetot[m2[1]+(1:m2[2]),m2[1]+(1:m2[2])] <- orthonormalization(base.exo)
    
    # Initialization
    # X associated projector computation
    ProjX <- interpol.matrix(m,m2) %*% basetot
    
    # Matrix of the linear form
    mat.rho <- diag(cst1/((1:m2[1])^2)+(1-cst1)/exp(1:m2[1]))
    mat.rho[1:nrow(d.rho),1:ncol(d.rho)] <- d.rho
    mat.rho.exo <- diag(cst1/((1:m2[2])^2)+(1-cst1)/exp(1:m2[2]))
    mat.rho.exo[1:nrow(d.rho.exo),1:ncol(d.rho.exo)] <- d.rho.exo
    rho <- matrix(0,nrow=sum(m2),ncol=sum(m2))
    rho[1:m2[1],1:m2[1]] <- mat.rho
    rho[1:nrow(d.a),m2[1]+(1:ncol(d.a))] <- d.a
    rho[m2[1]+(1:m2[2]),m2[1]+(1:m2[2])] <- mat.rho.exo
    
    # Coefficients
    Alpha <- matrix(0,sum(m2),sum(m2))
    Alpha[1:m2[1],1:m2[1]] <- diag(cst1/(1:m2[1]))
    Alpha[m2[1]+(1:m2[2]),m2[1]+(1:m2[2])] <- diag(cst1/(1:m2[2]))
    Alpha[1:nrow(alpha),1:ncol(alpha)] <- m2[1]*alpha
    Alpha[m2[1]+(1:ncol(alpha.conj)),1:nrow(alpha.conj)] <- m2[2]*t(alpha.conj)
    Alpha[1:nrow(alpha.conj),m2[1]+(1:ncol(alpha.conj))] <- m2[2]*alpha.conj

    # calcul de alpha.exo
    nmax <- max(dim(d.rho),dim(d.rho.exo),dim(alpha.conj),dim(d.a))
    mat1.d.rho <- matrix(0,nmax,nmax)
    mat1.d.rho.exo <- matrix(0,nmax,nmax)
    mat1.alpha.conj <- matrix(0,nmax,nmax)
    mat1.d.a <- matrix(0,nmax,nmax)
    mat1.d.rho[1:nrow(d.rho),1:ncol(d.rho)] <- d.rho
    mat1.d.rho.exo[1:nrow(d.rho.exo),1:ncol(d.rho.exo)] <- d.rho.exo
    mat1.d.a[1:nrow(d.a),1:ncol(d.a)] <- d.a
    mat1.alpha.conj[1:nrow(alpha.conj),1:ncol(alpha.conj)] <- alpha.conj
    alpha.exo <- invgen(mat1.d.rho.exo) %*% (mat1.alpha.conj -
                 mat1.d.rho.exo %*% mat1.alpha.conj %*% t(mat1.d.rho)) %*%
                     invgen(t(mat1.d.a))
    nrow.alpha.exo <- max((1:nrow(alpha.exo))[apply(alpha.exo!=0,1,any)])
    ncol.alpha.exo <- max((1:ncol(alpha.exo))[apply(alpha.exo!=0,2,any)])
    alpha.exo <- alpha.exo[1:nrow.alpha.exo,1:ncol.alpha.exo,drop=FALSE]

    # fin du calcul des coef   
    Alpha[m2[1]+(1:nrow(alpha.exo)),m2[1]+(1:ncol(alpha.exo))] <- m2[2]*alpha.exo
    Tmpmat <- Alpha - (rho %*% Alpha %*% t(rho))
    Tmpmat <- Tmpmat[-c((nrow(alpha)+1):m2[1],(nrow(alpha.exo)+1+m2[1]):sum(m2)),
                     ,drop=FALSE]
    Tmpmat <- Tmpmat[,-c((nrow(alpha)+1):m2[1],(nrow(alpha.exo)+1+m2[1]):sum(m2)),
                     drop=FALSE]
    test <- (try(Tmpmat <- t(pdMatrix(pdSymm(Tmpmat),factor=TRUE)),TRUE))
    if (inherits(test,"try-error"))
      stop("The resulting matrix of the noise is not positive definite.\n Modify d.rho, d.rho.exo and/or alpha, alpha.exo\n to get (alpha - (rho %*% alpha %*% t(rho))) positive definite.")
    Mu <- matrix(0,sum(m2),sum(m2))
    Mu[1:m2[1],1:m2[1]] <- diag(cst1/(1:m2[1]))
    Mu[m2[1]+(1:m2[2]),m2[1]+(1:m2[2])] <- diag(cst1/(1:m2[2]))
    Mu[1:nrow(alpha),1:nrow(alpha)] <- Tmpmat[1:nrow(alpha),1:nrow(alpha)]
    Mu[m2[1]+(1:nrow(alpha.exo)),m2[1]+(1:nrow(alpha.exo))] <-
        Tmpmat[nrow(alpha)+(1:nrow(alpha.exo)),nrow(alpha)+(1:nrow(alpha.exo))]
    Mu[1:nrow(alpha),m2[1]+(1:nrow(alpha.exo))] <-
        Tmpmat[1:nrow(alpha),nrow(alpha)+(1:nrow(alpha.exo))]
    Mu[m2[1]+(1:nrow(alpha.exo)),1:nrow(alpha)] <-
        Tmpmat[nrow(alpha)+(1:nrow(alpha.exo)),1:nrow(alpha)]
    
    # Strong white noise simulation
    epsilon <- matrix(rnorm(sum(m2)*n),nrow=sum(m2),ncol=n)
    epsilon <- Mu %*% epsilon
    
    # Creation of an empty X
    X <- matrix(0,nrow=sum(m2),ncol=n)
    
    # Iteration of simulation
    X[,1] <- epsilon[,1]
    for (i in (2:n))
        X[,i] <- rho %*% X[,i-1] + epsilon[,i]
    
    # Projection in the canonical basis
    restmp <- ProjX %*% X
    res <- list("X"=restmp[1:m[1],],"Z"=restmp[m[1]+1:m[2],])
    dimnames(res$X) <- list(paste(1:m[1]),paste(1:n))
    dimnames(res$Z) <- list(paste(1:m[2]),paste(1:n))
    return(as.fdata(res,name=c("X","Z")))
}

#  *****************************************************************************
#   Title : theoritical.coef
#         ************************************************************
#   Description : 
#       Calculation of the theoritical coefficient of a FARX model
#   Version : 0.8
#   Date : 2005-01-11
#  *****************************************************************************
theoretical.coef <- function(m=12,base=base.simul.far(24,5),
                       base.exo=NULL,
                       d.rho=diag(c(0.45,0.90,0.34,0.45)),
                       d.a=NULL,
                       d.rho.exo=NULL,
                       alpha=diag(c(0.5,0.23,0.018)),
                       alpha.conj=NULL,
                       cst1=0.05)
{
    # Testing the type (FAR or FARX)
    if (is.null(base.exo) || is.null(d.a)
        || is.null(d.rho.exo) || is.null(alpha.conj)) .farxmodel <- FALSE
    else .farxmodel <- TRUE

    if (.farxmodel) {
      if (length(m)==1) m <- rep(m,2)
      m2 <- c(nrow(base),nrow(base.exo))
    } else {
      m <- c(m,0)
      m2 <- c(nrow(base),0)
    }
  
    # Creation of the global basis
    if (.farxmodel) {
      basetot <- matrix(0,ncol=sum(m2),nrow=sum(m2))
      basetot[1:m2[1],1:m2[1]] <- orthonormalization(base)
      basetot[m2[1]+(1:m2[2]),m2[1]+(1:m2[2])] <- orthonormalization(base.exo)
    } else {
      basetot <- orthonormalization(base)
    }
    
    # Initialization
    # X associated projector computation
    ProjX <- interpol.matrix(m,m2) %*% basetot
    # Eigen vectors of X
    Xvect <- ProjX[1:m[1],1:ncol(base),drop=FALSE] 
    if (.farxmodel) {
      # Eigen vectors of Z
      Zvect <- ProjX[m[1]+(1:m[2]),m2[1]+(1:ncol(base.exo)),drop=FALSE]
    }
    
    # Matrix of the linear form
    mat.rho <- diag(cst1/((1:m2[1])^2)+(1-cst1)/exp(1:m2[1]))
    mat.rho[1:nrow(d.rho),1:ncol(d.rho)] <- d.rho
    rho <- matrix(0,nrow=sum(m2),ncol=sum(m2))
    rho[1:m2[1],1:m2[1]] <- mat.rho
    if (.farxmodel) {
      mat.rho.exo <- diag(cst1/((1:m2[2])^2)+(1-cst1)/exp(1:m2[2]))
      mat.rho.exo[1:nrow(d.rho.exo),1:ncol(d.rho.exo)] <- d.rho.exo
      rho[1:nrow(d.a),m2[1]+(1:ncol(d.a))] <- d.a
      rho[m2[1]+1:m2[2],m2[1]+1:m2[2]] <- mat.rho.exo
    }
    
    # Coefficients
    Alpha <- matrix(0,sum(m2),sum(m2))
    Alpha[1:m2[1],1:m2[1]] <- diag(cst1/(1:m2[1]))
    Alpha[1:nrow(alpha),1:ncol(alpha)] <- m2[1]*alpha
    if (.farxmodel) {
      Alpha[m2[1]+(1:m2[2]),m2[1]+(1:m2[2])] <- diag(cst1/(1:m2[2]))
      Alpha[m2[1]+(1:ncol(alpha.conj)),1:nrow(alpha.conj)] <- m2[2]*t(alpha.conj)
      Alpha[1:nrow(alpha.conj),m2[1]+(1:ncol(alpha.conj))] <- m2[2]*alpha.conj

      # calcul de alpha.exo
      nmax <- max(dim(d.rho),dim(d.rho.exo),dim(alpha.conj),dim(d.a))
      mat1.d.rho <- matrix(0,nmax,nmax)
      mat1.d.rho.exo <- matrix(0,nmax,nmax)
      mat1.alpha.conj <- matrix(0,nmax,nmax)
      mat1.d.a <- matrix(0,nmax,nmax)
      mat1.d.rho[1:nrow(d.rho),1:ncol(d.rho)] <- d.rho
      mat1.d.rho.exo[1:nrow(d.rho.exo),1:ncol(d.rho.exo)] <- d.rho.exo
      mat1.d.a[1:nrow(d.a),1:ncol(d.a)] <- d.a
      mat1.alpha.conj[1:nrow(alpha.conj),1:ncol(alpha.conj)] <- alpha.conj
      alpha.exo <- invgen(mat1.d.rho.exo) %*% (mat1.alpha.conj -
                   mat1.d.rho.exo %*% mat1.alpha.conj %*% t(mat1.d.rho)) %*%
                   invgen(t(mat1.d.a))
      nrow.alpha.exo <- max((1:nrow(alpha.exo))[apply(alpha.exo!=0,1,any)])
      ncol.alpha.exo <- max((1:ncol(alpha.exo))[apply(alpha.exo!=0,2,any)])
      alpha.exo <- alpha.exo[1:nrow.alpha.exo,1:ncol.alpha.exo,drop=FALSE]

      # end of the calculation
      Alpha[m2[1]+(1:nrow(alpha.exo)),m2[1]+(1:ncol(alpha.exo))] <- m2[2]*alpha.exo
    }

    dimX <- min(dim(alpha))
    valpX <- (eigen(alpha)$values)[1:dimX]
    if (.farxmodel) {
      dimZ <- min(dim(alpha.exo))
      dimT <- dimX + dimZ
      valpZ <- (eigen(alpha.exo[1:dimZ,1:dimZ])$values)[1:dimZ]
    } else dimT <- dimX
    eigenT <- eigen(Alpha)
    chgt.base <- eigenT$vectors[, 1:dimT,drop=FALSE]
    chgt.base2 <- matrix(0,ncol=dimT,nrow=sum(m2))
    chgt.base2[1:dimX,1:dimX] <- diag(dimX)
    if (.farxmodel) chgt.base2[m2[1]+(1:dimZ),dimX+(1:dimZ)] <- diag(dimZ)

    if (.farxmodel) {
      return(list("vectp.X"=Xvect[,1:dimX,drop=FALSE],
                  "vectp.Z"=Zvect[,1:dimZ,drop=FALSE],
                  "dim.X"=dimX, "dim.Z"=dimZ,
                  "valp.Var.X"=round(valpX,3), "valp.Var.Z"=round(valpZ,3),
                  "rho.X.Z"=round(t(chgt.base2)%*%rho%*%chgt.base2,3),
                  "vectp.T"=interpol.matrix(m,m2) %*% basetot %*% chgt.base,
                  "dim.T"=dimT,
                  "valp.Var.T"=round((eigenT$values)[1:dimT]/sum(m2),3),
                  "rho.T"=round(t(chgt.base)%*%rho%*%chgt.base,3)))
    } else {
      return(list("vectp.X"=Xvect[,1:dimX,drop=FALSE],
                  "dim.X"=dimX,
                  "valp.Var.X"=round(valpX,3),
                  "rho.X"=round(t(chgt.base2)%*%rho%*%chgt.base2,3)))
    }
}
