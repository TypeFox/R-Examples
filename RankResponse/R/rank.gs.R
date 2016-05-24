#    <rank.gs>
#    Copyright (C) <2014>  <Hsiuying Wang, Yu-Chun Lin>
#
#
#    This program is free software; you can redistribute it and/or modify
#
#    it under the terms of the GNU General Public License as published by
#
#    the Free Software Foundation; either version 2 of the License, or
#
#    (at your option) any later version.
#
# 
#
#    This program is distributed in the hope that it will be useful,
#
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#
#    GNU General Public License for more details.

rank.gs=function(data,alpha,type=2)
{
  data=as.matrix(data)
  data=data[!apply(apply(data,1,is.na),2,any),]
  n=dim(data)[1]
  k=dim(data)[2]
  z=qnorm(1-alpha/2)
  sumI=matrix(0,ncol=k,nrow=k)

  for(i in 1:(k-1))
  {
   a=sum(data[,i])/n
   for(j in (i+1):k)
   {
    b=sum(data[,j])/n
    if(type==1)
    {
       x=abs(a-b)/sqrt((a+b)/n)
    }else{

       c=sum(data[data[,i]==1&data[,j]==1,1])/n
       x=(a-b)/sqrt((a+b-2*c)/n)   
    }  
    if(x>=z)
    {
       sumI[i,j]=1
    }
   }
  }

  for(i in 2:k)
  {
   a=sum(data[,i])/n
   for(j in 1:(i-1))
   {
     b=sum(data[,j])/n
     if(type==1)
      {
         x=abs(a-b)/sqrt((a*(1-a)+b*(1-b)+2*a*b)/n)
      }else{

         c=sum(data[data[,i]==1&data[,j]==1,1])/n
         x=(a-b)/sqrt(((a-c)*(1-a+2*b-c)+(b-c)*(1-b+c))/n)   
      }  
      if(x>=z)
      {
        sumI[i,j]=1
      }
    }
  }


  rank=k-apply(sumI,1,sum)

  probability=apply(data,2,mean)
  result=rbind(probability,rank)

  return(result)

}
