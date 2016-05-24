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
### MVN
### iwish is in MCMCpack

RgetthetaMVN=function(y,z, priors1, priors2, theta, mus)  
{  
    K=dim(theta)[3]   #J,J,K
    J=dim(theta)[1]
   
    priors1=matrix(priors1, K,1)   
    theta=array(0,dim=c(J,J,K))
    for(k in 1:K)
    {   these=which(z==k)
        vv=length(these)
        sumy=matrix(0,J,J)
        for(tt in these)
        {   then=y[tt,]-mus[tt,] # w  ATJ
            sumy=sumy+then%*%t(then)
        }
        theta[,,k]=riwish(vv+priors1[k],sumy+priors2[,,k])   #2J df and 40*diag(J) are the IW prior for stability
    }
     
      theta
}
