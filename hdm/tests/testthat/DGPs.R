DGP.IV <-
function(n=250, p=100, pnz=50, Fstat=180, control=list(s2e=1, Cev=.6, s2z=1, szz=.5, pi1=1, alpha=1)) {
  s2e <- control$s2e
  Cev <- control$Cev
  s2z <- control$s2z
  szz <- control$szz
  pi1 <- control$pi1
  alpha <- control$alpha
  indp <- 0:(p-1)
  
  SZ <- s2z*toeplitz((szz)^indp)
  cSZ <- chol(SZ)
  cFS <- (c(rep(1,pnz), rep(0,p-pnz))*pi1)^indp
  scale <- matrix(0, nrow=1, ncol=1)
  s2v <- matrix(0, nrow=1, ncol=1)
  scale<- sqrt(Fstat/((Fstat+n)*t(cFS)%*%SZ%*%cFS))
  s2v <- 1-(scale^2)*t(cFS)%*%SZ%*%cFS
  sev <- Cev*sqrt(s2e)*sqrt(s2v)
  SU <- matrix(c(s2e,sev,sev,s2v), ncol=2)
  cSU <- chol(SU)
  
  zorig <- matrix(rnorm(n*p), nrow=n,ncol=p)%*%cSZ
  U <- matrix(rnorm(n*2), ncol=2)%*%cSU
  xorig <- scale[1,1]*zorig%*%cFS+U[,2]
  yorig <- alpha*xorig+U[,1]
  Z <- zorig - colMeans(zorig) #rep(1,n)*mean(zorig)
  X <- xorig - mean(xorig)
  Y <- yorig - mean(yorig)
  return(list(X=X, y=Y, Z=Z, setup=list(alpha=alpha, Pi=cFS, scale=scale, s2v=s2v, sev=sev, SU=SU)))
}

DGP.HC <- function(n=250, p=100, alpha=0.5, design=1, R2=c(0.5,0.5)) {
  Sigma <- toeplitz(0.5^(0:(p-1)))
  C <- chol(Sigma)
  X <- matrix(rnorm(n*p), nrow=n,ncol=p)%*%C
  beta01 <- 1/(1:p)^2
  beta02 <- 1/(1:p)^2 
  
  if (design==1) {
     sigma_d <- 1
     sigma_y <- 1
   }
   if (design==2) {
     sigma_d <- sqrt((1+X%*%beta01)^2/mean((1+X%*%beta01)^2))
   }
  
  cy <- sqrt(R2[1]/((1-R2[1])*t(as.matrix(alpha*beta02+ beta01))%*%Sigma%*%as.matrix(alpha*beta02 +  beta01)))
  beta01  <- cy * beta01
  cd <- sqrt(R2[2]/((1-R2[2])*t(as.matrix(beta02))%*%Sigma%*%as.matrix(beta02)))
  beta02  <- cd * beta02
  d <- X%*%beta02 + sigma_d*rnorm(n) 
  
  if (design==2) {
    sigma_y <-  sqrt((1+alpha*d+X%*%beta01)^2/mean((1+ alpha*d + X%*%beta01)^2))
  }


   y <- alpha*d+ X%*%beta01+ sigma_y*rnorm(n)

   return(list(y=y,X=X,d=d))
}

DGP.HCHIV <- function(n=250, px=300, pz=150, alpha=1) {
  Sigma <- toeplitz(0.5^(0:(px-1)))
  Sx <- chol(Sigma)
  Se <- chol(matrix(c(1, .6, .6,1), ncol=2))
  
  beta <- gamma <-  1/(1:px)^2
  delta <- 1/(1:pz)^2
  theta <- rbind(diag(pz), matrix(0, ncol=px-pz, nrow=pz))
  
  X <- matrix(rnorm(n*px), nrow=n, ncol=px)%*%Sx
  Z <- X%*%theta + 0.5*matrix(rnorm(n*pz), ncol=pz)
  
  e <- matrix(rnorm(n*2), ncol=2)
  
  d <- Z%*%delta + X%*%gamma + e[,1]
  y <- d*alpha + X%*%beta + e[,2]
  
  return(list(y=y, d=d, X=X, Z=Z))
}