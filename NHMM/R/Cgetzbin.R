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

Cgetzbin=function(K,z)
{   T=length(z)
    zbin=matrix(0,T,K)  #reset to zero
    for(k in 1:K)
    {  zbin[z==k,k]=1
    }
    zbin
}

getWbin=function(z,K,J)
{  T=length(z)
   
   Wbin1=array(0,dim=c(K,T,J))  #reset to zero
   for(k in 1:K)
   {  Wbin1[k,z==k,]=1
   }
   Wbin1
}