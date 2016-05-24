#-------------------------------------------------------------------------------
# November 2008
# Functions for cubic penalized regression spline which is represented as in
# equation (4.2) in Wood (2006): value-second-derivative (or cardinal) representation
#
# The functions are written based on the thoretical results shown in
# notes'CubicSplineNotes.pdf'
# (or Chapter 3 of the thesis...For now) and some of the functions are borrowed from 
# Wood 2006, Chapter 3.
#-------------------------------------------------------------------------------
# Bmat(), Dmat(), Xmat(), Smat(), trio.Xmat(), trio.Xmat.mul(), trio.Smat(), trio.Smat.mul() 
# Matrices Bmat and Dmat defines the relationship between
# beta_k (=f(x_k)) and delta_k (=f"(x_k)).
# The entries of Bmat and Dmat are as defined in Wood, pg 150 (2006).
Bmat <- function(xK){
  #xK is the vector of K knots
  # calculating h[k]=x[k+1]-x[k]
  K <- length(xK)
  h <- xK[-1] - xK[-K] #(h[1],h[2],...h[K])'
  Bmat <- matrix(0, nrow=(K-2), ncol=(K-2))
  for(i in 1:(K-2))
  diag(Bmat)[i] <- (h[i]+h[i+1])*(1/3)

  if(K>3){
    for(i in 1:(K-3)){
      Bmat[(i+1),i] <- Bmat[i,(i+1)] <- h[(i+1)]/6
    }
  }
  Bmat
}

Dmat <- function(xK){
  K <- length(xK)
  h <- xK[-1] - xK[-K]
  Dmat <- matrix(0, nrow=(K-2), ncol=K)
  for(i in 1:(K-2)){
    diag(Dmat)[i] <- 1/h[i]
    Dmat[i,(i+1)] <- -(1/h[i]+1/h[(i+1)])
    Dmat[i,(i+2)] <- 1/h[(i+1)]
  }
  Dmat
}

Xmat <- function(xK,x){
  #calculate the row vector of basis function at x.
  ##cat("range of x:")
  ##print(range(x))
  ##cat("range of xK:")
  ##print(range(xK))
  
  xK = as.vector(xK)
  n = length(x)
  K <- length(xK)
  if(any(xK[1]>x)|any(xK[K]<x))
    stop("some covariate values are not in the range of knots.")

  ##if(K>3){
  h <- xK[-1] - xK[-K]
  
  j.vec=NA
  ind = x==xK[1]
  j.vec[ind] = 1
  rm(ind)
  ind = x==xK[length(xK)]
  j.vec[ind] = (length(xK)-1)
  rm(ind)
  
  for(j in 2:length(xK)){
    ind = (x >= xK[(j-1)]) & (x < xK[j])
    j.vec[ind] = (j-1)
  }
  
  ##j.vec <- find.j2(x=x, xK=xK)
  
  ax.neg <- (xK[(j.vec+1)]-x)/h[j.vec]
  ax.pos <- (x-xK[j.vec])/h[j.vec]
  cx.neg <- ((xK[(j.vec+1)]-x)^3/h[j.vec] - h[j.vec]*(xK[(j.vec+1)]-x))/6
  cx.pos <- ((x-xK[j.vec])^3/h[j.vec] - h[j.vec]*(x-xK[j.vec]))/6
  
  a.vec <- cf.neg.vec <- cf.pos.vec <- matrix(0,ncol=K,nrow=n)
  
  ## calculating Fmat
  cBmat <- chol(Bmat(xK))
  Binv <- chol2inv(cBmat)
  Fmat <-  matrix(0, nrow=K, ncol=K)
  Fmat[2:(K-1),] <- Binv%*%Dmat(xK)
  ## end calculating Fmat
  
  for(j in 1:(K-1)){
    ind <- j.vec==j
    a.vec[ind,j] = ax.neg[ind]
    a.vec[ind,(j+1)] = ax.pos[ind]
    cf.neg.vec[ind,] <- kronecker(matrix(Fmat[j,],nrow=1),cx.neg[ind])
    cf.pos.vec[ind,] <- kronecker(matrix(Fmat[(j+1),],nrow=1),cx.pos[ind])
    ##cf.neg.vec[ind,] <- kronecker(matrix(Fmat(xK)[j,],nrow=1),cx.neg[ind])
    ##cf.pos.vec[ind,] <- kronecker(matrix(Fmat(xK)[(j+1),],nrow=1),cx.pos[ind])
    ##rm(ind)
  }
  
  X <- a.vec+cf.neg.vec+cf.pos.vec
  X
}

Smat <- function(xK){
  #penalty matrix 
  K <- length(xK)
  #cBmat <- chol(Bmat(xK)) #Bmat is positive-definite!!
  #Binv <- chol2inv(cBmat)
  #Smat <- t(Dmat(xK))%*%Binv%*%Dmat(xK) ## easy to calculate??
  Smat <- t(Dmat(xK))%*%solve(Bmat(xK))%*%Dmat(xK) ## easy to calculate??
  Smat
}

trio.Xmat <- function(x, mt, k1, k2, xk1, xk2, cenv){
  #x and mt can be from formatted data set 'triodat' obtained from pre.triogam()
  #non-genetic covariate x.
  #y1=1 if the risk allele had been transmitted to the affected child and 0, otherwise.
  #For mating type 3, (y1,y2); (0,0) if cgeno=++, (1,0) if +-, and (1,1) if --.
  #The basis is set up assuming that the basis for the first smoother is the same as
  #that of the second smoother; i.e., the lenghts and positions of the knots for the
  #two smoother are assumed to be the same.
  #rightnow, it is a little bit too slow;
  #should I use gam()'s smooth.construct.cr.smooth.spec() instead?
  n <-length(x) #i.e., n<-n+n(MT3)
  n1=length(unique(x[mt==1|mt==3.1]));##print(n1)
  n2=length(unique(x[mt==2|mt==3.2]));##print(n2)
  ind1 <- mt==1; ind2 <- mt==2; ind31 <- mt==3.1; ind32 <- mt==3.2

  if( !is.null(xk1) & !is.null(xk2) ){
    if ( (n1<k1)|(n2<k2) ) {
      if(n1<k1 & n2>=k2){ # not enough for GRR1
        msg <- paste(cenv, " has insufficient unique values to support ", 
                     k1, " knots: reduce k1.", sep = "")
        stop(msg)
      }
      
      else if(n1>=k1 & n2<k2){ # not enough for GRR2
        msg <- paste(cenv, " has insufficient unique values to support ", 
                     k2, " knots: reduce k2.", sep = "")
        stop(msg)
      }
      
      else{ # not enough for both GRR1 and GRR2
        msg <- paste(cenv, " has insufficient unique values to support ", 
                     k1, " and ", k2, " knots: reduce k1 and k2.", sep = "")
        stop(msg)
      }
      
    }
    
    ##gamma=1 #needed for estimating the smoothing parameter from GCV
    n.col <- k1 +k2 # for augmented X
    
    ##Constructing Design Matrix (i.e., cubic spline basis function matrix)
    
    ##Mating type 1
    X.tilde <- matrix(0, nrow=n, ncol=n.col)
    
    if(sum(ind1)==0) # 
      X.tilde <- X.tilde
    else 
      X.tilde[ind1, (1:k1)] <- Xmat(xK=xk1, x=x[ind1])

    ##Mating type 2
    if(sum(ind2)==0) # 
      X.tilde <- X.tilde
    else 
      X.tilde[ind2,(1+k1):(k1+k2)] <- Xmat(xK=xk2, x=x[ind2])

    ##Mating type 3
    if(sum(ind31)==0) # sum(ind31) must be equal to sum(ind32)
      X.tilde <- X.tilde
    else{ #x[ind31]=x[ind32]
      X.tilde[ind31, 1:k1] <- Xmat(xK=xk1, x=x[ind31])
      X.tilde[ind32, (1+k1):(k1+k2)] <- Xmat(xK=xk2, x=x[ind32])
    }#else ends
  }#if(!is.null(xk1) & !is.null(xk2))

  else {#either xk1 or xk2 is null
    if(!is.null(xk1)) {
      if(n1<k1){ # not enough for GRR1
        msg <- paste(cenv, " has insufficient unique values to support ", 
                     k1, " knots: reduce k1.", sep = "")
        stop(msg)
      }
      
      ##gamma=1 #needed for estimating the smoothing parameter from GCV
      ncol <- k1
      ##Constructing Design Matrix (i.e., cubic spline basis function matrix)
      
      ##Mating type 1
      X.tilde <- matrix(0, nrow=n, ncol=ncol)

      if(sum(ind1)==0) # 
        X.tilde <- X.tilde
      else
        X.tilde[mt==1,] <- Xmat(xK=xk1, x=x[ind1])
      
      ##Mating type 3
      if(sum(ind31)==0) # sum(ind31) must be equal to sum(ind32)
        X.tilde <- X.tilde
      else{ #x[ind31]=x[ind32]
        X.tilde[mt==3.1,] <- Xmat(xK=xk1, x=x[ind31])
      }#else ends
      
    } #if(is.null(xk1)) 

    else if(!is.null(xk2)) {
      if(n2<k2){ # not enough for GRR2
        msg <- paste(cenv, " has insufficient unique values to support ", 
                     k2, " knots: reduce k2.", sep = "")
        stop(msg)
      }
      ##gamma=1 #needed for estimating the smoothing parameter from GCV

      ncol <- k2
      ##Constructing Design Matrix (i.e., cubic spline basis function matrix)
      
      ##Mating type 2
      X.tilde <- matrix(0, nrow=n, ncol=ncol)
      
      if(sum(ind2)==0) # 
        X.tilde <- X.tilde
      else 
        X.tilde[mt==2,] <- Xmat(xK=xk2, x=x[ind2])

      ##Mating type 3
      if(sum(ind32)==0) # sum(ind31) must be equal to sum(ind32)
        X.tilde <- X.tilde
      else{ #x[ind31]=x[ind32]
        X.tilde[mt==3.2,] <- Xmat(xK=xk2, x=x[ind32])
      }
    } #else if(is.null(xk2))
    
    else { #both xk1 and xk2 are null
      stop("knots are null for both GxE functions")
    } #else: both xk1 and xk2 are null
  } # else: either xk1 or xk2 is null
  
  X.tilde
}

trio.Xmat.mul <- function(x, xk, cenv){
  #x and mt can be from formatted data set 'triodat' obtained from pre.triogam()
  #non-genetic covariate x.
  #y1=1 if the risk allele had been transmitted to the affected child and 0, otherwise.
  #For mating type 3, (y1,y2); (0,0) if cgeno=++, (1,0) if +-, and (1,1) if --.
  #The basis is set up assuming that the basis for the first smoother is the same as
  #that of the second smoother; i.e., the lenghts and positions of the knots for the
  #two smoother are assumed to be the same.
  #rightnow, it is a little bit too slow;
  #should I use gam()'s smooth.construct.cr.smooth.spec() instead?
  n <-length(x); k <- length(xk) #i.e., n<-n+n(MT3)
  stopifnot(n>k)
  n.col <- k 
  
  ##Constructing Design Matrix (i.e., cubic spline basis function matrix)
  
  ##Mating type 1
  X.tilde <- matrix(0, nrow=n, ncol=n.col)
  X.tilde <- Xmat(xK=xk, x=x)
  X.tilde
}

trio.Smat <- function(xk1,xk2,pen.scale,X,mt){#sp1,sp2,
 # have to make sure xk1 and xk2 are the same as those used for 
 # constructing X1 and X2!!

  ind1 <- mt==1 | mt==3.1
  ind2 <- mt==2 | mt==3.2

  k1 <- length(xk1); k2 <- length(xk2)
  
  
  if(k1>0 & k2>0){#both xk1 and xk2 are non null
    X1 <- X[ind1,1:k1]
    X2 <- X[ind2,(k1+1):(k1+k2)]
    pen.mat1 <- Smat(xk1)
    pen.mat2 <- Smat(xk2)
  }#both xk1 and xk2 are non null
  else if(k1>0 & k2==0){
    X1 <- X; X2 <- NULL
    pen.mat1 <- Smat(xk1); pen.mat2 <- NULL
  }# xk1 is not null
  else if(k1==0 & k2>0){
    X1 <- NULL; X2 <- X
    pen.mat <- NULL; pen.mat2 <- Smat(xk2)
  } #else if(k1==0 & k2>0
  else {#both k1 and k2 are zero
    stop("knots are null for both GxE functions")
  } #both k1 and k2 are bigger than zero 

  dim.smat <- k1+k2
  pen.mat <- matrix(0, nrow=dim.smat, ncol=dim.smat)
  
  if(pen.scale){
    ## taken from smoothCon(): 
    ## "The following is intended to make scaling `nice' for better gamm performance.
    ##  Note that this takes place before any resetting of the model matrix, and 
    ##  any `by' variable handling. From a `gamm' perspective this is not ideal, 
    ##  but to do otherwise would mess up the meaning of smoothing parameters
    ##  sufficiently that linking terms via `id's would not work properly (they 
    ##  would have the same basis, but different penalties)
    ##  if (scale.penalty && length(sm$S)>0 && is.null(sm$no.rescale)) # then the penalty coefficient matrix is rescaled
    ## 
    ##  {  maXX <- mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X
    ##      for (i in 1:length(sm$S)) {
    ##        maS <- mean(abs(sm$S[[i]]))
    ##        sm$S[[i]] <- sm$S[[i]] * maXX / maS
    ##      }  
    ##  }"
    if(!is.null(X1)){
      maXX <- mean(abs(t(X1) %*% X1))
      maS <- mean(abs(pen.mat1))
      pen.scale1 <- (maXX/maS)
      pen.mat1 <- pen.mat1*pen.scale1
    }

    if(!is.null(X2)){
      maXX <- mean(abs(t(X2) %*% X2))
      maS <- mean(abs(pen.mat2))
      pen.scale2 <- (maXX/maS)
      pen.mat2 <- pen.mat2*pen.scale2
    }
  }#if(pen.scale)
  
  else{
    pen.mat1 <- pen.mat1
    pen.mat2 <- pen.mat2
  }
  
  ## multiply pen matrices by sp
  ##pen.mat1 <- sp1*pen.mat1
  ##pen.mat2 <- sp2*pen.mat2

  if( (k1>0) & (k2>0) ){
    pen.mat[1:k1,1:k1] <- pen.mat1
    pen.mat[(k1+1):(k1+k2),(k1+1):(k1+k2)] <- pen.mat2
  }

  else if( (k1>0) & (k2==0) ){
    pen.mat <- pen.mat1
  }

  else{
    pen.mat <- pen.mat2
  }
  
  S <- pen.mat
  
  S
}
trio.Smat.mul <- function(xk,X,pen.scale=TRUE){#sp1,sp2,
 # have to make sure xk1 and xk2 are the same as those used for 
 # constructing X1 and X2!!

  k <- length(xk)
  stopifnot(k>0)
  pen.mat <- Smat(xk)
  
  if(pen.scale){
    ## taken from smoothCon(): 
    ## "The following is intended to make scaling `nice' for better gamm performance.
    ##  Note that this takes place before any resetting of the model matrix, and 
    ##  any `by' variable handling. From a `gamm' perspective this is not ideal, 
    ##  but to do otherwise would mess up the meaning of smoothing parameters
    ##  sufficiently that linking terms via `id's would not work properly (they 
    ##  would have the same basis, but different penalties)
    ##  if (scale.penalty && length(sm$S)>0 && is.null(sm$no.rescale)) # then the penalty coefficient matrix is rescaled
    ## 
    ##  {  maXX <- mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X
    ##      for (i in 1:length(sm$S)) {
    ##        maS <- mean(abs(sm$S[[i]]))
    ##        sm$S[[i]] <- sm$S[[i]] * maXX / maS
    ##      }  
    ##  }"
    if(!is.null(X)){
      maXX <- mean(abs(t(X) %*% X))
      maS <- mean(abs(pen.mat))
      pen.scale <- (maXX/maS)
      pen.mat <- pen.mat*pen.scale
    }
  }#if(pen.scale)
  
  else{
    pen.mat <- pen.mat
  }
  pen.mat  
}
