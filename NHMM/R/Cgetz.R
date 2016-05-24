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


Cgetz=function( z, QQ, denzity, subseqy)   
{   T=dim(denzity)[2]
    K=dim(denzity)[1]
    pro=numeric(K)
    qs=matrix(0,K,T)
    
    t=1
    qs[,1]=QQ[,z[1+1],1+1]  #t=1
    denzity1=matrix(0,K,T)
    for(k in 1:K)
    {  denzity1[k,]=apply(denzity[k,,],1,prod)
    }  
       
    pro=denzity1[,t]*qs[,t]                #probabilities
    if(sum(pro)==0)  #TRUE - find any zero probs
    {  pro=rep(1/K,K)               #if all probabilities are zero, set all states as equally likely                   
    }
    z[t]=rcpp_rmultinom(pro) 
    
    
    for(t in 2:(T-1))
    {      if(subseqy[t-1] != subseqy[t]  )   #check for sequence breaks
           {  qs[,t]=QQ[,z[t+1],t+1]
           }else{  if(subseqy[t] != subseqy[t+1])
                   {  qs[,t]=QQ[z[t-1],,t]
                   }else{  qs[,t]=QQ[z[t-1],,t]*QQ[,z[t+1],t+1]  #no special case (majority of the time)
                   }	
           } 				
   
    
        pro=denzity1[,t]*qs[,t]                #probabilities
        if(sum(pro)==0)  #TRUE - find any zero probs
        {  pro=rep(1/K,K)               #if all probabilities are zero, set all states as equally likely
        }
        z[t]=rcpp_rmultinom(pro)          #z[t]=sample( 1:K, 1 ,prob=pro)
        
    }
    
    t=T
    qs[,T]=QQ[,z[T-1],T]  #T=1
    pro=denzity1[,t]*qs[,t]                #probabilities
    if(sum(pro)==0)  #TRUE - find any zero probs
    {  pro=rep(1/K,K)               #if all probabilities are zero, set all states as equally likely
    }
    z[t]=rcpp_rmultinom(pro) 

    z
}