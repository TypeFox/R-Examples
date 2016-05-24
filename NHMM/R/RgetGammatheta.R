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

### Gamma (nmix, W) not used
## delta =FALSE
## no W
## no nmix



RgetGammatheta=function(y,z, priors, theta, nmix, vvv,delt)  
{  
  K=dim(theta)[3]   #, nmix, K,J
  J=dim(theta)[4]
  
  if(is.na(priors[2,1,1,1]))  #fixed AA parameter
  { for(v in 1:nmix)
    {AA=matrix(theta[1,v,,],K,J)   #K by J
     BB=matrix(theta[2,v,,],K,J)
     a.AA=matrix(priors[1,v,,],K,J)   #4,nmix,K,J
     b.AA=matrix(priors[2,v,,],K,J)    #precision
     a.BB=matrix(priors[3,v,,],K,J)
     b.BB=matrix(priors[4,v,,],K,J)   #K by J
     prop=matrix(priors[5,v,,],K,J)
     
     for(k in 1:K)
     {  for(j in 1:J)
        {  n=sum(z==k & vvv[,j]==(v-1+delt))
           if(n > 0) #ensure there is data in this state, if not skip
           {  yy=(y[,j])[z==k & vvv[,j]==(v-1+delt)]
              theta[1,v,k,j]=a.AA[k,j]
              theta[2,v,k,j]=rgamma(1,n*theta[1,v,k,j] + a.BB[k,j], b.BB[k,j]+sum(yy))
          
             #print("uhoh")
             #print(yy)
           }
        }
     }
  }
  theta
  }else{
    for(v in 1:nmix)
    {  AA=matrix(theta[1,v,,],K,J)   #K by J
       BB=matrix(theta[2,v,,],K,J)
       a.AA=matrix(priors[1,v,,],K,J)   #4,nmix,K,J
       b.AA=matrix(priors[2,v,,],K,J)    #precision
       a.BB=matrix(priors[3,v,,],K,J)
       b.BB=matrix(priors[4,v,,],K,J)   #K by J
       prop=matrix(priors[5,v,,],K,J)
       
       lll=numeric(5)
       for(k in 1:K)
       {  for(j in 1:J)
          {   n=sum(z==k & vvv[,j]==(v-1+delt))
              if(n > 1) #ensure there is data in this state, if not skip
              {  yy=(y[,j])[z==k & vvv[,j]==(v-1+delt)]
                 theta[2,v,k,j]=rgamma(1,n*AA[k,j] + a.BB[k,j], b.BB[k,j]+sum(yy))
                 aprop=rnorm(1,AA[k,j],prop[k,j]) #propopsal
                 if(aprop > 0)
                 {  #one=sum(log(dgamma(yy,aprop,b.AA[k,j])*dgamma(aprop,a.AA[k,j], b.AA[k,j])))
                    #two=sum(log(dgamma(yy,AA[k,j],b.AA[k,j])*dgamma(AA[k,j],a.AA[k,j], b.AA[k,j])))
                    one=sum(log(dgamma(yy,aprop,b.AA[k,j])/dgamma(yy,AA[k,j],b.AA[k,j])))
                    two=sum(log(dgamma(aprop,a.AA[k,j], b.AA[k,j])/dgamma(AA[k,j],a.AA[k,j], b.AA[k,j])))
                    
                    alpha.mh=exp(one+two)
                   #one=(log(dgamma(yy,aprop,b.AA[k,j])))
                   #two=(log(dgamma(yy,AA[k,j],b.AA[k,j])))
                   #three=(log(dgamma(aprop,a.AA[k,j], b.AA[k,j])))
                   # four=(log(dgamma(AA[k,j],a.AA[k,j], b.AA[k,j])))
                   #alpha.mh=exp(sum(one-two+three-four))

                 
                    if(runif(1) < alpha.mh){ theta[1,v,k,j]=aprop }
                 }
              }
           }
       }
    }
    theta
  }
}



