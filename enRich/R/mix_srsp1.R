mix_srsp1 <-function(datasr, datasp, method=NULL, initialpara_sr, initialpara_sp, stopdiff=0.0001, Nsr, Nsp)
{
	if (length(datasr)>0)
	{
		N=nrow(datasr)
	}else
	{
		N=nrow(datasp)
	}

	if (length(Nsr)==1)
	{
		if(Nsr==0)
		Ns=0
		else
		Ns=length(Nsr)
	}else
	{
		Ns=length(Nsr)
	}
	Np=Ns+Nsp
	temp1=0
	para_sr=NULL
	para_sp=NULL
	for (i in 1: sum(Nsr))
	{
		if (length(initialpara_sr)>0)
		{
			if (method=="Poisson")
			{
				temp1<-temp1+initialpara_sr[(i-1)*4+1]
				para_sr=c(para_sr, initialpara_sr[((i-1)*4+2):((i-1)*4+4)])
			}
			if (method=="NB")
			{
				temp1<-temp1+initialpara_sr[(i-1)*6+1]
				para_sr=c(para_sr, initialpara_sr[((i-1)*6+2):((i-1)*6+6)])
			}
		}
	}
	for (i in 1: sum(Nsp))
	{
		if (length(initialpara_sp)>0)
		{
			if (method=="Poisson")
			{
				temp1<-temp1+initialpara_sp[(i-1)*4+1]
				para_sp=c(para_sp, initialpara_sp[((i-1)*4+2):((i-1)*4+4)])
			}
			if (method=="NB")
			{
				temp1<-temp1+initialpara_sp[(i-1)*6+1]
				para_sp=c(para_sp, initialpara_sp[((i-1)*6+2):((i-1)*6+6)])
			}
		}
	}	
   	p=temp1/(sum(Nsr)+sum(Nsp))
   	para=c(p, para_sr, para_sp)
   	k1=1
   	stopN=10000
   	difference=rep(1, length(para))

	if (Nsp==1)
	datasp=cbind(datasp, 0)

   	#Iterative steps 
   	while (any(difference>stopdiff)&(k1<stopN))
   	{
      		paratemp=para
  	 
		#1 iterative classification (E-step)
	 	pr1=matrix(0, N, Np)
   		pr0=matrix(0, N, Np)
	   	ppr1=matrix(0, N, Np)
   		ppr0=matrix(0, N, Np)

		if (Ns>0)
		{
			for (j in 1:Ns)## same regions for technical replicates
			{
             		for (j1 in 1:Nsr[j])
				{
                   			j2=sum(Nsr[0:(j-1)])+j1
					if (method=="Poisson")
					{	
	            				pr1[,j]=pr1[,j]+ifelse(datasr[,j2]<para_sr[(j2-1)*3+3], log(0.00001), dpois(datasr[,j2]-para_sr[(j2-1)*3+3],para_sr[(j2-1)*3+1], log=TRUE))
			   			pr0[,j]=pr0[,j]+dpois(datasr[,j2], para_sr[(j2-1)*3+2], log=TRUE)
					}
					if (method=="NB")
					{	
	      		      			pr1[,j]=pr1[,j]+ifelse(datasr[,j2]<para_sr[(j2-1)*5+5], log(0.00001), dnbinom(datasr[,j2]-para_sr[(j2-1)*5+5], para_sr[(j2-1)*5+2],,para_sr[(j2-1)*5+1], log=TRUE))
				   		pr0[,j]=pr0[,j]+dnbinom(datasr[,j2], para_sr[(j2-1)*5+4],,para_sr[(j2-1)*5+3], log=TRUE)
					}
             	 	}
			}
		}
		if (Nsp>0)
		{
			for (j in (Ns+1):Np)## same p for non-replicates
			{	
				j2=j-Ns
				if (method=="Poisson")
				{
					pr1[, j]=ifelse(datasp[,j2]<para_sp[(j2-1)*3+3], log(0.00001), dpois(datasp[,j2]-para_sp[(j2-1)*3+3], para_sp[(j2-1)*3+1], log=TRUE))
	   				pr0[, j]=dpois(datasp[,j2], para_sp[(j2-1)*3+2], log=TRUE)
				}
				if (method=="NB")
				{
					pr1[, j]=ifelse(datasp[,j2]<para_sp[(j2-1)*5+5], log(0.00001), dnbinom(datasp[,j2]-para_sp[(j2-1)*5+5], para_sp[(j2-1)*5+2],,para_sp[(j2-1)*5+1], log=TRUE))
	   				pr0[, j]=dnbinom(datasp[,j2], para_sp[(j2-1)*5+4],,para_sp[(j2-1)*5+3], log=TRUE)
				}
			}
		}        
		for (j in 1:Np)
		{
             	pr1[,j]=pr1[,j]+log(p)
			pr0[,j]=pr0[,j]+log(1-p)
   			ppr1[,j]=1/(1+exp(pr0[,j]-pr1[,j]))
	      		ppr0[,j]=1-ppr1[,j]
		}
      		#2. Maximum likelihood estimates (M-step)
	      	# optim/optimize function
      		p=sum(ppr1)/(N*Np)
		if (Ns>0)
		{
		      	for (j in 1:Ns)
			{
				for (j1 in 1:Nsr[j])
				{
 					j2=sum(Nsr[0:(j-1)])+j1
					if (method=="Poisson")
					{
				      		para_sr[(j2-1)*3+1]=sum(ppr1[,j]*datasr[,j2])/sum(ppr1[,j])-para_sr[(j2-1)*3+3] #signal mean
					      	para_sr[(j2-1)*3+2]=sum(ppr0[,j]*datasr[,j2])/sum(ppr0[,j]) #background mean
					}
					if (method=="NB")
					{
				      		para_sr[(j2-1)*5+1]=sum(ppr1[,j]*datasr[,j2])/sum(ppr1[,j])-para_sr[(j2-1)*5+5] #signal mean
					      	para_sr[(j2-1)*5+3]=sum(ppr0[,j]*datasr[,j2])/sum(ppr0[,j]) #background mean
						temp=optim(c(para_sr[(j2-1)*5+2], para_sr[(j2-1)*5+4]), fn=NBlike_phi_offset, data=datasr[,j2], ppr1=ppr1[,j], ppr0=ppr0[,j], mus=para_sr[(j2-1)*5+1], mub=para_sr[(j2-1)*5+3], k=para_sr[(j2-1)*5+5])$par
 		     				para_sr[(j2-1)*5+2]=temp[1]
		   	   			para_sr[(j2-1)*5+4]=temp[2]
					}
				}
			}
		}
		if (Nsp>0)
		{
			for (j in (Ns+1):Np)## same p for non-replicates
			{	
				j2=j-Ns
				if (method=="Poisson")
				{
					para_sp[(j2-1)*3+1]=sum(ppr1[,j]*datasp[,j2])/sum(ppr1[,j])-para_sp[(j2-1)*3+3] #signal mean
				      	para_sp[(j2-1)*3+2]=sum(ppr0[,j]*datasp[,j2])/sum(ppr0[,j]) #background mean
				}
				if (method=="NB")
				{
					para_sp[(j2-1)*5+1]=sum(ppr1[,j]*datasp[,j2])/sum(ppr1[,j])-para_sp[(j2-1)*5+5] #signal mean
			      		para_sp[(j2-1)*5+3]=sum(ppr0[,j]*datasp[,j2])/sum(ppr0[,j]) #background mean
					temp=optim(c(para_sp[(j2-1)*5+2], para_sp[(j2-1)*5+4]), fn=NBlike_phi_offset, data=datasp[,j2], ppr1=ppr1[,j], ppr0=ppr0[,j], mus=para_sp[(j2-1)*5+1], mub=para_sp[(j2-1)*5+3], k=para_sp[(j2-1)*5+5])$par
		 	     		para_sp[(j2-1)*5+2]=temp[1]
   					para_sp[(j2-1)*5+4]=temp[2]
				}
			}
		}
		para=c(p, para_sr, para_sp)
	      difference=abs(para-paratemp)
      		k1=k1+1
   	}

	## subfunction
	arrpara<-function(para, method)
	{
		if (method=="Poisson")
		{
			Nexp=length(para)/3
			parameter=matrix(0, Nexp, 3)
			colnames(parameter)=c("lambda_S", "lambda_B", "k")
			for (j in 1:Nexp)
	   		{
				parameter[j,1:3]=para[((j-1)*3+1):((j-1)*3+3)]
			}
	   	}
		if (method=="NB")
		{
			Nexp=length(para)/5
	   		parameter=matrix(0, Nexp, 5)
	   		colnames(parameter)=c("mus", "phis", "mub", "phib", "k")
		   	for (j in 1:Nexp)
   			{
				parameter[j,1:5]=para[((j-1)*5+1):((j-1)*5+5)]
	   		}
		}
	   	rownames(parameter)=rownames(parameter,  do.NULL = FALSE, prefix = "Experiment")
   		return(parameter)
	}
	if (Ns>0)
	{
		temp_sr=arrpara(para_sr, method=method)
	   	temp_sr=cbind(rep(p,sum(Nsr)), temp_sr)
		colnames(temp_sr)[1]="p" 	
		temp_sr=round(temp_sr, digits=4)
	}else
	{
		temp_sr=NULL
	}
	if (Np>0)
	{
		temp_sp=arrpara(para_sp, method=method)
   		temp_sp=cbind(rep(p,sum(Nsp)), temp_sp)
		colnames(temp_sp)[1]="p"
		temp_sp=round(temp_sp, digits=4)
	}else
	{
		temp_sr=NULL
	}
   	parameter=list(para_sr=temp_sr, para_sp=temp_sp)
	return(parameter)
}
