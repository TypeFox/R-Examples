#-------------------------------------------------------
#
#  Created       : 30/10/02
#  Last Modified : Time-stamp: <2003-04-02 09:52:47 lucas>
#
#  Description   : Robust principal component analysis
#                  
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#
#  Licence       : GPL 
#
#-------------------------------------------------------

K <- function(u,kernel="gaussien") {
    switch(kernel,
        gaussien = (2*pi)^(-1/2) * exp(-u^2/2),
        quartic   = 15/16 * (1-u^2)^2 * (abs(u)<1),
        triweight = 35/32 * (1-u^2)^3 * (abs(u)<1),
        epanechikov = 3/4 * (1-u^2) *   (abs(u)<1),
        cosinus = pi/4 * cos (u*pi/2) * (abs(u)<1),
        uniform = 1/2 * (abs(u)<1),
    )
}

# Variance locale
W <- function(x,h,D=NULL,kernel="gaussien")
{
    x   <- as.matrix(x)
    n   <- dim(x)[1]
    p   <- dim(x)[2]
    if (is.null(D)) {
        D <- diag(1,p)
    }
    x <- as.vector(x)
    D <- as.vector(D)
    kernel <- substr(kernel,1,1)

    VarLoc <- .C(
                 "W",
                 as.double(x),
                 as.double(h),
                 as.double(D),
                 as.integer(n),
                 as.integer(p),
                 as.character(kernel),
                 res=double(p*p),
                 result = as.integer(1),
                 PACKAGE= "amap"
                 )

    if(VarLoc$result == 2)
      stop("Cannot allocate memory")
    if(VarLoc$result == 1)
      stop("Error")

    matrix(VarLoc$res,p)
}


varrob <- function(x,h,D=NULL,kernel="gaussien")
{
    x   <- as.matrix(x)
    x   <- scale(x, center = TRUE, scale = FALSE)
    n   <- dim(x)[1]
    p   <- dim(x)[2]
    if (is.null(D)) {
        D <- diag(1,p)
    }
    x <- as.vector(x)
    D <- as.vector(D)
    kernel <- substr(kernel,1,1)

    Calcul <- .C(
                 "VarRob",
                 as.double(x),
                 as.double(h),
                 as.double(D),
                 as.integer(n),
                 as.integer(p),
                 as.character(kernel),
                 res=double(p*p),
                 result = as.integer(1),
                 PACKAGE= "amap")

    if(Calcul$result == 2)
      stop("Cannot allocate memory")
    if(Calcul$result == 1)
      stop("Error")

    S <- matrix(Calcul$res,p)
    Sinv <- solve(S)
    solve ( Sinv - D / h)
}


acpgen <- function(x,h1,h2,center=TRUE,reduce=TRUE,kernel="gaussien")
{
    # CENTRONS, ET REDUISONS
    x    <- as.matrix(x)
    x    <- scale(x ,center = center, scale = FALSE)
    if (reduce == TRUE)
         {
          x    <- apply(x,2,function(u) { u/sd(u)}) 
         }

    # ESTIMATION DE W et VarRob
    n <- dim(x)[1]
    VarInv   <- solve(var(x)*(n-1)/n) # solve= inverser
    leU    <- varrob(x,h1,D=VarInv,kernel=kernel)
    leW    <- W(x,h2,D=VarInv,kernel=kernel)
    Winv   <- solve(leW) 


    # anal. spec de Var.W^-1 :
    EIG    <- eigen(leU %*% Winv)  
    V      <- EIG$vector

    #EIG    <- eigen( x %*% Winv %*% t(x)  )
    #U      <- EIG$vector
    #n      <- dim(x)[1]
    #p      <- dim(x)[2]
    #S      <- diag(Re(EIG$values),n)   
    #S1     <- diag(Re(1/EIG$values),n)
    #S      <- sqrt(S[,1:p])
    #S1     <- sqrt(S1[,1:p])
    #V      <- t(x)%*% U%*% S1
    # X=U.S.V' -> V = X' U S^-1
    

    # AFFICHAGE DES RESULTATS


    scores <- x %*% Winv %*% V

    V      <- as.matrix(V)
    scores <- as.matrix(scores)
    dimnames(V)[[2]] <- paste("Comp",1:dim(x)[2])
    if(!is.null( dimnames(x)[[2]] ))
      dimnames(V)[[1]] <- dimnames(x)[[2]]
    if(!is.null( dimnames(x)[[1]] ))
      dimnames(scores)[[1]] <- dimnames(x)[[1]]
    dimnames(scores)[[2]] <- paste("Comp",1:dim(x)[2])
    eig    <- sqrt(EIG$values)
    sdev   <- apply(scores,2,sd)    
    res    <- list(eig=eig,sdev=sdev,scores=scores,loadings=V)
    class(res) <- "acp"
    res
}


acprob <- function(x,h=1,center=TRUE,reduce=TRUE,kernel="gaussien")
{   
    x    <- as.matrix(x)
    x    <- scale(x ,center = center, scale = FALSE)
    if (reduce == TRUE)
         {
          x    <- apply(x,2,function(u) { u/sd(u)}) 
         }
    EIG  <- eigen( varrob(x,h),symmetric=TRUE) 
    V    <- EIG$vector    # ou bien: V=svd(x)$v

    val  <- sqrt(EIG$values)

    scores <- x %*% V

    V      <- as.matrix(V)
    scores <- as.matrix(scores)
    dimnames(V)[[2]] <- paste("Comp",1:dim(x)[2])
    if(!is.null( dimnames(x)[[2]] ))
      dimnames(V)[[1]] <- dimnames(x)[[2]]
    if(!is.null( dimnames(x)[[1]] ))
      dimnames(scores)[[1]] <- dimnames(x)[[1]]
    dimnames(scores)[[2]] <- paste("Comp",1:dim(x)[2])
    sdev   <- apply(scores,2,sd)    
    res  <- list(eig=val,sdev=sdev,scores=scores,loadings=V)
    class(res) <- "acp"
    res
}
