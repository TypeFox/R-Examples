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


RgetPoissontheta=function(y,z, priors, theta, nmix, vvv,delt)  
{  
  K=dim(theta)[3]   #, nmix, K,J
  J=dim(theta)[4]
  
  for(v in 1:nmix)
  {   for(v in 1:nmix)
      {  
         a.AA=matrix(priors[1,v,,], K,J)   #4,nmix,K,J
         b.AA=matrix(priors[2,v,,], K,J)    #precision
         #as.matrix(priors[3,v,,])
         #as.matrix(priors[4,v,,])   #K by J
         #as.matrix(priors[5,v,,])
         
         for(k in 1:K)
         {  for(j in 1:J)
            {  n=sum(z==k & vvv[,j]==(v-1+delt))
               if(n > 0) #ensure there is data in this state, if not skip
               {      theta[1,v,k,j]=rgamma(1, a.AA[k,j]+sum((y[,j])[z==k & vvv[,j]==(v-1+delt)]), b.AA[k,j]+n)
                      theta[2,v,k,j]=NA
               }
            }
          }
      }
  }
  theta
}