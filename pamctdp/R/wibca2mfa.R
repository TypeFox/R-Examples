#----------------------------------------------------------------------------------
# Calculo de los dos AFMTC asociados al ACIBP 
#----------------------------------------------------------------------------------
# Funcion wibca2mfa 
# Agosto 30 de 2010
# Campo Elias Pardo
# SALE
#------------------------------------------------------------------------------
# ==========funcion para calcular Lg
   Lg <- function(Xj,Mj,D=diag(nrow(Xj))/nrow(Xj),Xk=Xj,Mk=Mj)
   {
     Wj <- Xj %*% Mj   %*% t( Xj)
     Wk <- Xk %*% Mk   %*% t( Xk)
     Lg <- sum(diag(Wj%*%D %*%Wk %*%D))
     return(Lg)
   }  
#=======================================
#========== funcion ibca2mfa
wibca2mfa <- function(ACww)
{
  rbl <- ACww$rbl
  cbl <- ACww$cbl
  J <- nrow(ACww$cbvar)
  L <- nrow(ACww$lbvar)
  nf <- ACww$nf
  # valores propios
  # afm considerando bandas columna
  homJ <- ACww$hom[2]
  colb <- NULL
  colb$eig<-ACww$eig*homJ
  # afm considerando bandas columnafila
  homL <- ACww$hom[1]
  rowb <- NULL
  rowb$eig<-ACww$eig*homL
  # coordenadas de las bandas
  colb$coor  <-homJ*(ACww$cbvar * (ACww$cbw %*% t(rep(1,nf))))
  rowb$coor  <-homL*(ACww$lbvar * (ACww$lbw %*% t(rep(1,nf))))
  # Lg
   X <- as.matrix(ACww$tab)
   cbl.fac <- rep(1:J,cbl)
   rbl.fac <- rep(1:L,rbl)
   D <- diag(ACww$lw)
   M<- diag(ACww$cw)
   # bandas columna
   Mbc <- M*homJ
   ev1 <- eigen(t(X)%*%D%*%X%*%M,symmetric=FALSE, only.values = TRUE)$values[1]
   Lgbc <- matrix(NA,J+1,J+1)
   for (j in 1:J)
     {
       Xj <- X[,cbl.fac==j];Mj <- Mbc[cbl.fac==j,cbl.fac==j]
       Lgbc[J+1,j]<- Lg(Xj,Mj,D,X,M/ev1)
       for (k in 1:J)
         {
           Xk <- X[,cbl.fac==k];Mk <- Mbc[cbl.fac==k,cbl.fac==k]
           Lgbc[j,k]  <- Lg(Xj,Mj,D,Xk,Mk)
         }
     }
   Lgbc[J+1,J+1]<- Lg(X,M/ev1,D)
   Lgbc[,J+1]<- Lgbc[J+1,]
   rownames(Lgbc)<-colnames(Lgbc)<- c(rownames(ACww$cbvar),"MFA")
   RVbc <- (diag(1/(sqrt(diag(Lgbc)))))%*%Lgbc%*% (diag(1/sqrt(diag(Lgbc))))
   rownames(RVbc)<-rownames(Lgbc); colnames(RVbc)<-colnames(Lgbc)
  
   #----------sale------------
   colb$Lg <- Lgbc
   colb$RV <- RVbc
   #--------------------------         
   # bandas fila
   #--------------------------
   Dbl <- D*homL 
   Lgbl <- matrix(NA,L+1,L+1)
   for (l in 1:L)
            {
              Xl <- t(X[rbl.fac==l,])
              Dl <- Dbl[rbl.fac==l,rbl.fac==l]
              Lgbl[L+1,l]<- Lg(Xl,Dl,M,t(X),D/ev1)
              for (m in 1:L)
                {
                  Xm <- t(X[rbl.fac==m,]);Dm <- Dbl[rbl.fac==m,rbl.fac==m]
                  Lgbl[l,m]  <- Lg(Xl,Dl,M,Xm,Dm)
                }
            }
   Lgbl[L+1,L+1]<- Lg(t(X),D/ev1,M)
   Lgbl[,L+1]<- Lgbl[L+1,]
   rownames(Lgbl)<-colnames(Lgbl)<- c(rownames(ACww$lbvar),"MFA")
   RVbl <- (diag(1/sqrt(diag(Lgbl))))%*%Lgbl%*% (diag(1/sqrt(diag(Lgbl))))
   rownames(RVbl)<-rownames(Lgbl); colnames(RVbl)<-colnames(Lgbl)
   #----------sale------------
   rowb$Lg <- Lgbl
   rowb$RV <- RVbl
   #--------------------------
   mfa <- list(colb=colb,rowb=rowb)       
   return(mfa)
 } # fin de la funcion
 #================================================= 
