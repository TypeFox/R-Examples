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

#library(msm) rtnorm
################################# STARTS HERE
RgetM=function(A,K, psi, gamy, Wbin, gamboo, vvv )
{  T=dim(vvv)[1]
   J=dim(vvv)[2]
   M=matrix(0,T,J)   
   mmu=matrix(0,T,J)
   LL=dim(psi)[1]
   
   for(j in 1:J)
   {  mmu=matrix(psi[,j],1,LL) %*% Wbin[,,j] #beta0[z+1,j]+ENSO1m[,j]*beta1[j]+ENSO2m[,j]*beta3[j]
      #M[,j]=Crtnorm(T, mmu, rep(1.0,T), gams[vvv[,j]+1], gams[vvv[,j]+2], gamboo[vvv[,j]+1], gamboo[vvv[,j]+2] )        
       M[,j]=rtnorm(T,mmu, rep(1,T), gamy[vvv[,j]+1,j], gamy[vvv[,j]+2,j])
   }
   M
}

   