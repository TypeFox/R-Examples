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
#############################################################
Rgetpsi=function(M, psi, Wbin, psipriorm, psipriorp, A, K)
{   J=dim(M)[2]
    T=dim(M)[1]
    
    Wj=matrix(0,K+A,T)
    PL.j=matrix(0,K+A,K+A)
    bL.j=matrix(0,K+A,1)
    P1.j=matrix(0,K+A,K+A)
    V1.j=matrix(0,K+A,K+A)
    m1.j=matrix(0,K+A,1)
    sqrtV1.j=matrix(0,K+A,K+A)
    

    for(j in 1:J)
    {  psij=matrix(psi[,j],1, A+K)  #A+K
       Wj=matrix(Wbin[,,j],A+K,T)   #A+K by T
  
       if(!is.element(0,apply(Wj,1,sum)))  ### are all K are used?
       {   PL.j = Wj %*%  t(Wj);  
           bL.j = Wj %*% matrix(M[,j],T,1)
           P1.j = PL.j + diag(A+K)*psipriorp[,j];  #add the variance prior
           V1.j = chol2inv(chol(P1.j));     
           m1.j = V1.j %*% (bL.j) + psipriorp[,j]*psipriorm[,j];   #add the mean prior
           sqrtV1.j = t(chol(V1.j));
           psi[,j] = m1.j + sqrtV1.j %*% rnorm(A+K);
       }else{
          a=which(apply(Wj,1,sum)==0)
          K2=K-length(a)
          Wj2=Wj[-a,]
          psi2=as.matrix(psi[-a,])
          psipriorm2=psipriorm[-a,]
          psipriorp2=psipriorp[-a,]
          
           PL.j = Wj2 %*%  t(Wj2);  
           bL.j = Wj2 %*% as.matrix(M[,j])
           P1.j = PL.j + diag(A+K2)*psipriorp2[,j];  #add the variance prior
           V1.j = chol2inv(chol(P1.j));     
           m1.j = V1.j %*% (bL.j) + psipriorp2[,j]*psipriorm2[,j];   #add the mean prior
           sqrtV1.j = t(chol(V1.j));
           psi[-a,j] = m1.j + sqrtV1.j %*% rnorm(A+K2);
       }
          
    }
    psi
}

