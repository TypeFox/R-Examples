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
RgetpsiMVN=function(y, z, psi, Wbin, psipriorm, psipriorp, A, K, thetainv)
{   J=dim(y)[2]
    T=dim(y)[1]
    sigsig=matrix(0,J,J)
    
    thetainv1=thetainv
    ##include the prior
    for(k in 1:K)
    {  thetainv1[,,k]=thetainv[,,k]+psipriorp[k,]*diag(J)
    }
    
    ####  psi0 terms  (first 1:K terms)  
    for(k in 1:K)
    {  these=which(z==k)
       if(length(these)>0)
       {  sumz1=matrix(0,J,J);  sumz2=numeric(J)
          for(tt in these)
          {   sumz1=sumz1+thetainv1[,,k]
              if(A==1){Y2=y[tt,]-psi[(K+1):(K+A),]*Wbin[(K+1):(K+A),tt,]}
              if(A>1){Y2=y[tt,]-apply(psi[(K+1):(K+A),]*Wbin[(K+1):(K+A),tt,],2,sum)}
              if(A==0){Y2=y[tt,]}
              sumz2=sumz2+  Y2%*% thetainv1[,,z[tt]]
          }
          sigsig=solve(sumz1)
          mumu=t(sigsig)%*%t(sumz2+psipriorm[k,]*psipriorp[k,])
          psi[k,]=mvrnorm(1,mumu, sigsig)
       }
    }
    
    if(A>0) ## update the psi for the input variables
    {  for(ll in 1:A)
       {  sumz1=rcpp_getsumz1(K,J,T,z-1, c(thetainv1), Wbin[K+ll,,])      
          sumz2=rcpp_getsumz2(ll,A,K,J,T,z-1, c(thetainv1), c(Wbin[(K+1):(K+A),,]), y, t(psi[(K+1):(K+A),]), matrix(psi[1:K,],K,J)) 
          sigsig=solve(sumz1)
          mumu=t(sigsig)%*%(sumz2+matrix(psipriorm[(K+1):(K+A),]*psipriorp[(K+1):(K+A),],A,J)[ll,])
          psi[K+ll,]=mvrnorm(1,mumu, sigsig)
       }
    }
    
    psi
}

