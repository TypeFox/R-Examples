################################################################
## Copyright 2014 Tracy Holsclaw.

## This file is part of NHMM.

## NHMM is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or any later version.

## NHMM is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## NHMM.  If not, see <http://www.gnu.org/licenses/>.

## This code may have similar features as BayesLogit package but is adapted to the NHMM setting.
#############################################################



Rgetbeta=function( zbin,beta,XX, betapriorm, betapriorp)
{   T=dim(zbin)[1]
    K=dim(zbin)[2]
    L=dim(XX)[2]
    wij=diag(T) 
    C=matrix(0,T,K)
    

  if(!is.element(0,apply(zbin,2,sum)))  ### all z's are used
  {
      for(k in 1:(K-1))  #1:(K-1) without any skipped z
      {  C[,k]=log(apply(exp(beta[-k,]%*% t(XX)),2,sum))  #hold one row of beta out at a time
         wij=diag(T)*rpg(T,1,t(beta[k,]%*%t(XX))-C[,k])   #library
         
         PL.j = t(XX) %*% wij %*% XX;   #37 secs per 1000
         bL.j = t(XX) %*% (zbin[,k]-0.5+wij %*% C[,k]);
         P1.j = PL.j + betapriorp[,,k];  #add the variance prior
         V1.j = chol2inv(chol(P1.j));        
         m1.j = V1.j %*% (bL.j) + betapriorp[,,k]%*%betapriorm[,k];   #add the mean prior
         sqrtV1.j = t(chol(V1.j));
         beta[k,] = m1.j + sqrtV1.j %*% rnorm(L);  #beta is K by L
      } 
  }else{   #remove unused z components
    
     #warning("One of your states disappeared: think of setting K smaller")
     a=which(apply(zbin,2,sum)==0)
     K2=K-length(a)
     zbin2=zbin[,-a]
     beta2=as.matrix(beta[-a,-a])
     XX2=XX[,-a]
     betapriorm2=betapriorm[-a,-a]
     betapriorp2=betapriorp[-a,-a,-a]
     C2=matrix(0,T,K2)
     
     for(k in 1:(K2-1))  #1:(K-1) without any skipped z
     {  C2[,k]=log(apply(exp(beta2[-k,]%*% t(XX2)),2,sum))  #hold one row of beta out at a time
        wij=diag(T)*rpg(T,1,t(beta2[k,]%*%t(XX2))-C2[,k])   #library
        
        PL.j = t(XX2) %*% wij %*% XX2;   #37 secs per 1000
        bL.j = t(XX2) %*% (zbin2[,k]-0.5+wij %*% C2[,k]);
        P1.j = PL.j + betapriorp2[,,k];  #add the variance prior
        V1.j = chol2inv(chol(P1.j));        
        m1.j = V1.j %*% (bL.j) + betapriorp2[,,k]%*%betapriorm2[,k];   #add the mean prior
        sqrtV1.j = t(chol(V1.j));
        beta2[k,] = m1.j + sqrtV1.j %*% rnorm(L-length(a));  #beta is K by L
     } 
     beta[-a,-a]=beta2
  }

   beta
}