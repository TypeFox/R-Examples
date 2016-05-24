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
### Normal


RgetNormaltheta=function(y,z, priors, theta, nmix, vvv,delt)  
{  
    K=dim(theta)[3]   #, nmix, K,J
    J=dim(theta)[4]
  
      for(v in 1:nmix)
      {   a.mu=matrix(priors[1,v,,], K,J)   #4,nmix,K,J
          b.mu=matrix(priors[2,v,,], K,J)    #precision
          a.sig=matrix(priors[3,v,,], K,J)
          b.sig=matrix(priors[4,v,,], K,J)   #K by J
          #as.matrix(priors[5,v,,])
     
          for(k in 1:K)
          {   for(j in 1:J)
              {  n=sum(z==k & vvv[,j]==(v-1+delt))
                 if(n > 0) #ensure there is data in this state, if not skip
                 {  theta[2,v,k,j]=1/rgamma(1,n/2+a.sig[k,j], 1/2*sum((y[z==k & vvv[,j]==(v-1+delt),j]-theta[1,v,k,j])^2)+b.sig[k,j])
                    ss=1/(n/theta[2,v,k,j]+b.mu[k,j])
                    theta[1,v,k,j]=rnorm(1,(sum(y[z==k & vvv[,j]==(v-1+delt),j]/theta[2,v,k,j]+a.mu[k,j]*b.mu[k,j]))*ss,sd=sqrt(ss))
                 }
              }
           }
      }
      theta
}
