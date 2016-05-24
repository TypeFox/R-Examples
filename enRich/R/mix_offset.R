mix_offset <-function(data1, method=NULL, initialpara=NULL, k=0, stopdiff=0.0001)
{
## data1 is the vector of counts of the data
## method= "Poisson" or "NB"

## sufunction
#Log-likelihood based on mixture of NB(mu, size), maximised using optim() function
	NBlike_phi_offset<-function(para, data1, ppr1, ppr0, mus,mub, k=k)
	{
	   N=length(data1)
	   phis=para[1]
	   phib=para[2]
	   temp1=sum(ppr1*ifelse(data1<k, 0, dnbinom(data1-k, phis, ,mus, log=TRUE)))
	   temp2=sum(ppr0*dnbinom(data1, phib, ,mub, log=TRUE))
	   logl=temp1+temp2
	   return(-logl)
	}

	N=length(data1)
	para=initialpara
	k1=1
 	stopN=10000
   	difference=rep(1, length(para))
	  
	#2. Iterative steps 
   	while (any(difference>stopdiff)&(k1<stopN))
   	{
 	     	paratemp=para
	      #2.1 iterative classification (E-step)
		if (method=="Poisson")
		{
      			pr1=log(para[1])+ifelse(data1<k, log(0.00001), dpois(data1-k, para[2], log=TRUE))
		      	pr2=log(1-para[1])+dpois(data1, para[3], log=TRUE)
		}
		
		if (method=="NB")
		{
			pr1=log(para[1])+ifelse(data1<k, log(0.00001), dnbinom(data1-k, para[3], ,para[2], log=TRUE))
		      	pr2=log(1-para[1])+dnbinom(data1, para[5],,para[4], log=TRUE)    
		}
	 	ppr1=1/(1+exp(pr2-pr1))
		ppr0=1-ppr1

	      #2.2. Maximum likelihood estimates (M-step) 
		para[1]=sum(ppr1)/N
 
		if (method=="Poisson")
		{
			para[2]=sum(ppr1*data1)/sum(ppr1)-k
			para[3]=sum(ppr0*data1)/sum(ppr0)
		}	
		if (method=="NB")
		{
		     	para[2]=sum(ppr1*data1)/sum(ppr1)-k
		     	para[4]=sum(ppr0*data1)/sum(ppr0)	
		    	# optim/optimize function
      	     		temp=optim(c(para[3],para[5]), fn=NBlike_phi_offset, data=data1, ppr1=ppr1, ppr0=ppr0, mus=para[2], mub=para[4], k=k)$par
                 	para[3]=temp[1]
                 	para[5]=temp[2]
		}
	      difference=abs(para-paratemp)
		k1=k1+1
	}
	para=c(para, k)
	if (method=="Poisson")
   	names(para)=c("p", "lambda_S", "lambda_B", "offset")
      	if (method=="NB")
      	names(para)=c("p", "mu_S", "phi_S", "mu_B", "phi_B", "offset")
   	para=round(para, digits=4)
	return(para)

}
