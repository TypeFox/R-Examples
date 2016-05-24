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


Rgetgams=function(gamy, vvv,M,nmix,mixes)
{  
   if(mixes>2)
   { J=dim(vvv)[2]
     for(i in 1:(mixes-2))  #start at 2 because first gamy is set to zero
     {  for(j in 1:J)
        {  if(  sum(vvv[,j]==(i+1))!=0  & sum(vvv[,j]==i )!=0)   #dont update if no data
           {  gamy[i+2,j]=runif(1,max(0,max(M[vvv[,j]==i,j])),min(M[(vvv[,j])==(i+1),j]))
           }
        } 
     }
   }
  gamy
  
}   
