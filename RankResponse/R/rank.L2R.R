#    <rank.L2R>
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

rank.L2R=function(data,response.number,prior.parameter,e)
{
      data=as.matrix(data)
      v=response.number
      alpha=prior.parameter
      A <- matrix(0,2^v,v)     
      for(j in 1:v)
      {
		A[,j] <- rep(c(rep(0,2^(v-j)),rep(1,2^(v-j))),2^(j-1))
	}

      Temp=rbind(data,A)
      group=kmeans(Temp,(2^v))
      obs=numeric(2^v)
      for(j in 1:(2^v))
      {
           for(i in 1:(2^v))
           {
              if(all(group$center[j,]==A[i,])==T){obs[i]=(group$size[j])-1}
           }
      }

      n =(v-1)*v               
      a1 = numeric(n)   
      a2 = numeric(n)
      for(i in 1:v)
      {
		a1[((i-1)*(v-1)+1):((i-1)*(v-1)+(v-1))] = rep(i,(v-1))
	}

      a2[1:(v-1)] = c(2:v)
	a2[((v-1)*(v-1)+1):((v-1)*(v-1)+(v-1))] = c(1:(v-1))
	for(i in 2:(v-1)){
		a2[((v-1)*(i-1)+1):((v-1)*(i-1)+(v-1))] = c(1:(i-1),(i+1):v)
	}

	AA = cbind(A,alpha,obs)  

	V = matrix(0,1,n)       
		
	for(i in 1:n){

		alpha0 = sum(AA[,(v+1)]+AA[,(v+2)])       		
		alphai = AA[((AA[,a1[i]]!= AA[,a2[i]])& AA[,a1[i]]==1),][,(v+1)]+AA[((AA[,a1[i]]!=AA[,a2[i]])& AA[,a1[i]]==1),][,(v+2)]
		alphaj = AA[((AA[,a1[i]]!= AA[,a2[i]])& AA[,a2[i]]==1),][,(v+1)]+AA[((AA[,a1[i]]!=AA[,a2[i]])& AA[,a2[i]]==1),][,(v+2)]

		E = sum(alphai-alphaj)/alpha0             
		
		Var1 = sum(alphai*alpha0-alphai^2)
		Var2 = sum(alphaj*alpha0-alphaj^2)
		Var3 = 2* sum(alphai%*%t(alphaj))

		Var = (Var1+Var2+Var3)/(alpha0^3+alpha0)  
		V[1,i] = pnorm(-E/sqrt(Var),0,1)          
	}

	tv = round(V[1,],6)
	pt = cbind(a1,a2,tv)

      pi=numeric(v)
      N=sum(AA[,(v+2)])
      for(i in 1:v){

            pi[i]=sum(AA[AA[,i]==1,(v+2)])/N
      }
      
      names(pi)=c(1:v)
      pi_temp=sort(pi,decreasing=T)
      rank_temp=as.numeric(names(pi_temp))    
            
      u=numeric(v-1)
      for(i in 1:(v-1)){
            
            u[i]=pt[pt[,1]==rank_temp[i] & pt[,2]==rank_temp[i+1],3]
      }


      epsilon = 0.001		
	FD_bar = numeric(1) ; FDR_bar = numeric(100)

      for(j in 1:100){
		t = 0.01*j
		D = numeric(1) ; FD_bar = numeric(1)
		for(i in 1:(v-1)){
			if(u[i] > t){
				FD_bar = FD_bar + (1-u[i])
				D = D+1
			}
		}

		FDR_bar[j] = FD_bar/(D+epsilon)
	}

	for(i in 1:100){
		if(FDR_bar[i] <= e){break}
		t = i+1				
	}

      tD_star = t*0.01

      sumI_temp=!(tv>=tD_star)
      sumI=matrix(sumI_temp,ncol=v-1,nrow=v,byrow=T)
      rank=v-apply(sumI,1,sum)

      probability=apply(data,2,mean)
      result=rbind(probability,rank)

      return(result)

}