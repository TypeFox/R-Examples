
D_lm <- function(F, l, m, Beta)
{
  p <- dim(F)[1]
  lm <- max(l,m)
  q <- dim(F)[3]-1
  D <- matrix(0,p,p)
  for(i in 1:(p-1)){
   for(j in (i+1):p){
    for(k in (-q+lm):(q-lm)){
      D[i,j] <- D[i,j]+F[i,i,abs(k+l)+1]*F[j,j,abs(k+m)+1]+
                F[i,i,abs(k+l)+1]*F[j,j,abs(k-m)+1] 
    }
    D[i,j] <- 0.5*D[i,j]+0.25*(Beta[i,j]-1)*((F[,,l+1]+
           t(F[,,l+1]))[i,j]*(F[,,m+1]+t(F[,,m+1]))[i,j])              
   }
  }

  D <- D+t(D)
  
  for(i in 1:p){
   for(k in (-q+lm):(q-lm)){
    D[i,i] <- D[i,i]+F[i,i,abs(k+l)+1]*F[i,i,abs(k+m)+1]+
                     F[i,i,abs(k+l)+1]*F[i,i,abs(k-m)+1]
   }
   D[i,i] <- D[i,i]+(Beta[i,i]-3)*F[i,i,l+1]*F[i,i,m+1]
  }
  D
}


ASCOV_SOBI <- function(psi, taus, a=2, Beta=NULL, A=NULL)
{
  p <- dim(psi)[2]
  q <- dim(psi)[1]
  K <- length(taus)  

  if(is.null(A)) A <- diag(p)

  if(q<(3*K)) psi <- rbind(psi,matrix(0,3*K-q,p))
  
  q <- dim(psi)[1]

  Psi <- matrix(0,q+1,p)
  for(i in 1:p){
    Psi[,i] <- c(1,psi[,i])
    Psi[,i] <- Psi[,i]/sqrt(sum(Psi[,i]^2))
 }
  
  PSI <- array(0,c(p,p,q+1)) 
  for(i in 1:p){
   for(j in 1:(q+1)){
     PSI[i,i,j] <- Psi[j,i]
   }  
  }

  F_tau <- array(0,c(p,p,q+1))
  for(i in 0:q){
   for(j in 1:(q+1-i)){
      F_tau[,,i+1] <- F_tau[,,i+1]+tcrossprod(diag(PSI[,,j]),diag(PSI[,,j+i]))
   }
  }

  if(is.null(Beta)) Beta <- 2*diag(p)+matrix(1,p,p)

  Lambda <- array(0,c(p,p,K))
  for(k in 1:K){
   for(j in 1:p){  
    Lambda[j,j,k] <- F_tau[j,j,taus[k]+1]
   }
  }
    
  Sum_lam <- rep(0,p)
  for(j in 1:p){
    for(k in 1:K){
      Sum_lam[j] <- Sum_lam[j]+Lambda[j,j,k]^2
    }
  }   

  P <- matrix(0,p,p)
  ord <- order(Sum_lam,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]] <- 1
  }  
  
  for(k in 1:K){
    Lambda[,,k] <- tcrossprod(crossprod(t(P),Lambda[,,k]),P)
  }

  for(i in 1:(q+1)){
    F_tau[,,i] <- tcrossprod(crossprod(t(P),F_tau[,,i]),P)
  }

  ASCOV<-matrix(0,p^2,p^2)
  for(j in 1:p){
   for(i in 1:p){
    if(i!=j){
      ASV <- .C("ascov", as.double(as.vector(F_tau)),as.double(as.vector(Lambda)), as.double(taus),as.integer(c(i-1,j-1,p,q,K)),as.double(as.vector(Beta)),as.double(a),res=double(2), PACKAGE="BSSasymp")$res
 
      ASCOV <- ASCOV+ASV[2]*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(diag(p)[,j],diag(p)[,i]))+ASV[1]*kronecker(tcrossprod(diag(p)[,j],diag(p)[,j]),tcrossprod(diag(p)[,i],diag(p)[,i]))
    }  

    if(i==j) ASCOV <- ASCOV+0.25*D_lm(F_tau,0,0,Beta)[i,i]*kronecker(tcrossprod(diag(p)[,i],diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))  
   }
  }    

  EMD <- sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W <- crossprod(t(P),solve(A))
  W <- crossprod(diag(sign(rowMeans(W))),W)
  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)
}


# M is the number of autocovariances used in estimation.

ASCOV_SOBI_estN <- function(X, taus, mixed=TRUE, M=100, a=2)
{
  p <- dim(X)[2]
  T <- dim(X)[1]
  K <- length(taus)  
  
  if(mixed){  
    W <- SOBI(X,taus)$W
  }else W <- diag(p)

  X <- tcrossprod(sweep(X,2,colMeans(X)),W)    
  
  F_tau <- array(0,c(p,p,M+1))
  for(m in 0:M){
     F_tau[,,m+1] <- tcrossprod(t(X[1:(T-m),]),t(X[(m+1):T,]))/(T-m)
  }

  Beta <- 2*diag(p)+matrix(1,p,p)

  Lambda <- array(0,c(p,p,K))
  for(k in 1:K){
   for(j in 1:p){  
    Lambda[j,j,k] <- F_tau[j,j,taus[k]+1]
   }
  }

  Sum_lam <- rep(0,p)
  for(j in 1:p){
    for(k in 1:K){
      Sum_lam[j] <- Sum_lam[j]+Lambda[j,j,k]^2
    }
  }   

  P <- matrix(0,p,p)
  ord <- order(Sum_lam,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]] <- 1
  }  
  
  
  for(k in 1:K){
    Lambda[,,k] <- tcrossprod(crossprod(t(P),Lambda[,,k]),P)
  }

  for(i in 1:(M+1)){
    F_tau[,,i] <- tcrossprod(crossprod(t(P),F_tau[,,i]),P)
  }

  W <- crossprod(t(P),W)

  ASCOV <- matrix(0,p^2,p^2)
  for(j in 1:p){
   for(i in 1:p){
    if(i!=j){
       ASV <- .C("ascov", as.double(as.vector(F_tau)),as.double(as.vector(Lambda)), as.double(taus),as.integer(c(i-1,j-1,p,M,K)),as.double(as.vector(Beta)),as.double(a),res=double(2), PACKAGE="BSSasymp")$res
 
       ASCOV <- ASCOV+ASV[2]*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(diag(p)[,j],diag(p)[,i]))+ASV[1]*kronecker(tcrossprod(diag(p)[,j],diag(p)[,j]), tcrossprod(diag(p)[,i],diag(p)[,i]))
   }  

   if(i==j) ASCOV <- ASCOV+0.25*D_lm(F_tau,0,0,Beta)[i,i]*kronecker(tcrossprod(diag(p)[,i],diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))
  }
 }    

  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/T  
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/T
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A)
}


ASCOV_SOBI_est <- function(X, taus, arp=NULL, maq=NULL, mixed=TRUE, M=100, a=2, ...)
{
  p <- dim(X)[2]
  T <- dim(X)[1]
  K <- length(taus)  
 
  if(mixed){
    W <- SOBI(X,taus)$W
  }else W <- diag(p)
    
  X <- tcrossprod(sweep(X,2,colMeans(X)),W)  

  if(is.null(arp)) arp <- rep(1,p) 
  if(is.null(maq)) maq <- rep(1,p)

  Psi <- matrix(0,M+1,p)
  for(i in 1:p){
    arma <- arima(X[,i],c(arp[i],0,maq[i]),...)
    psi0 <- ARMAtoMA(ar=arma$coef[min(1,arp[i]):arp[i]],
              ma=arma$coef[(arp[i]+1):(arp[i]+maq[i])],M)
    Psi[,i] <- c(1,psi0)
    Psi[,i] <- Psi[,i]/sqrt(sum(Psi[,i]^2))
  }

  Beta <- matrix(0,p,p)
  for(i in 1:p){
    Ex4 <- mean(X[,i]^4)
    SumPsi4 <- sum(Psi[,i]^4)
    Beta[i,i] <- (Ex4-3)/SumPsi4+3
  }  
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Exi2xj2 <- mean(X[,i]^2*X[,j]^2)
      SumPsii2j2 <- sum(Psi[,i]^2*Psi[,j]^2)
      Beta[i,j] <- (Exi2xj2-1)/SumPsii2j2+1    
      Beta[j,i] <- Beta[i,j]
    }  
  }

  PSI <- array(0,c(p,p,M+1)) 
  for(i in 1:p){
   for(j in 1:(M+1)){
     PSI[i,i,j] <- Psi[j,i]
   }  
  }

  F_tau <- array(0,c(p,p,M+1))
  for(i in 0:M){
   for(j in 1:(M+1-i)){
      F_tau[,,i+1] <- F_tau[,,i+1]+tcrossprod(diag(PSI[,,j]),diag(PSI[,,j+i]))
   }
  }

  
  for(m in 0:M){
     diag(F_tau[,,m+1]) <- diag(tcrossprod(t(X[1:(T-m),]),t(X[(m+1):T,]))/(T-m))
  }
  
  Lambda <- array(0,c(p,p,K))
  for(k in 1:K){
   for(j in 1:p){  
    Lambda[j,j,k] <- F_tau[j,j,taus[k]+1]
   }
  }
  
  Sum_lam <- rep(0,p)
  for(j in 1:p){
    for(k in 1:K){
      Sum_lam[j] <- Sum_lam[j]+Lambda[j,j,k]^2
    }
  }   

  P <- matrix(0,p,p)
  ord <- order(Sum_lam,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]] <- 1
  }  
  
  for(k in 1:K){
    Lambda[,,k] <- tcrossprod(crossprod(t(P),Lambda[,,k]),P)
  }

  for(i in 1:(M+1)){
    F_tau[,,i] <- tcrossprod(crossprod(t(P),F_tau[,,i]),P)
  }

  W <- crossprod(t(P),W)

  ASCOV <- matrix(0,p^2,p^2)
  
  for(j in 1:p){
    for(i in 1:p){
      if(i!=j){
        ASV <- .C("ascov", as.double(as.vector(F_tau)),as.double(as.vector(Lambda)), as.double(taus),as.integer(c(i-1,j-1,p,M,K)),as.double(as.vector(Beta)),as.double(a),res=double(2), PACKAGE="BSSasymp")$res
 
        ASCOV <- ASCOV+ASV[2]*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(diag(p)[,j],diag(p)[,i]))+ASV[1]*kronecker(tcrossprod(diag(p)[,j],diag(p)[,j]),tcrossprod(diag(p)[,i],diag(p)[,i]))
      }  

      if(i==j) ASCOV <- ASCOV+0.25*D_lm(F_tau,0,0,Beta)[i,i]*kronecker(tcrossprod(diag(p)[,i],diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))   
    }
  }    

  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/T  
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/T
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A)
}



