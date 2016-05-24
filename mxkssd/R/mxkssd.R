#m=number of factors
#n=number of runs
#level_vec=a vector of order k containing the levels of factors such that number of factors having these levels are same.
#for example a design with 2^7,3^14 will have level_vec=c(2,3,3) so that 7 factors each having levels 2, 3 and 3 respectively.
#D is the supersaturated design

genvec=function(m,n,level_vec,k)
#The function created initial k-circulant q-level SSD randomly.
{
	for (i in 1:k)
	{
		if (n%%level_vec[i] !=0) stop("n should be divisible by levels of the factors")
	}		
	if (m!=(n-1)*k) stop("k-circulant design not possible") #checks condition 1	
	gvec=matrix(0,1,m)
	gmat=matrix(0,n-1,k)
	for (i1 in 1:k)
	{
		q=level_vec[i1]
		t=n/q
		for (i in 1:(q-1))
		{
			for (j in 1:t)
			{
				gmat[j+(i-1)*t,i1]=i+1
			}
		}
		for (i in (t*(q-1)+1):(n-1))
		{
			gmat[i,i1]=1
		}
		gmat[,i1]=sample(gmat[,i1],(n-1))
		#condition 2
	}		
	count=0
	for (j in 1:(n-1))
	{
		for (i in 1:k)
		{
			count=count+1
			gvec[,count]=gmat[j,i]
		}
	}	
	gvec
}

D_from_G=function(genvec,m,n,k)
{
#the function creates the design obtained by k-circulating the generator vector
	D=matrix(0,n,m)
	D[1,]=genvec
	#k-circulating the first row to next (n-2) rows
	for (rown in 1:(n-2))
	{
		for (coln in 1:m)
		{
			coln_plusk=coln+k
			if ((coln+k)>m) coln_plusk= (coln+k) %% m
			D[rown+1,coln_plusk]=D[rown,coln]
		}
	}	
	#Allocation of level q to last row of the design
	for (i in 1:m) D[n,i]=1
	D
}


fNOD=function(D,m,n,level_vec,k)
{
	#The function calculates average chisquare of the design	
	sum_fNOD=0	
	max_fNOD=0	
	for (i in 1:(m-1))
	{
		factor1=i%%k
		if (i%%k==0) factor1=k
		qi=level_vec[factor1]		
		for (j in (i+1):m)
		{
			factor2=j%%k
			if (j%%k==0) factor2=k
			qj=level_vec[factor2]
			freq_levels=matrix(-n/(qi*qj),qi,qj) #freq_levels is the matrix to contain freq of level pairs in the columns i and j
			for (i1 in 1:n)
			{
				freq_levels[D[i1,i],D[i1,j]]=freq_levels[D[i1,i],D[i1,j]]+1				
			}
			for (a in  1:qi)
			{
				for (b in 1:qj)
				{
					freq_levels[a,b]=(freq_levels[a,b])**2										
				}				
			}			
			temp1=sum(freq_levels)			
			sum_fNOD=sum_fNOD+temp1			
			if (temp1>max_fNOD) max_fNOD=temp1						
		}
	}	
	E_fNOD=sum_fNOD*2/(m*(m-1))	
	output1=list(E_fNOD,max_fNOD)
	output1
}

alias=function(D,m,n,level_vec,k)
{
	#The function calculates average chisquare of the design
	aliased=0		
	for (i in 1:(m-1))
	{
		factor1=i%%k
		if (i%%k==0) factor1=k
		qi=level_vec[factor1]
		freq_levels1=matrix(-n/(qi*qi),qi,qi)
		for (j in (i+1):m)
		{
			factor2=j%%k
			if (j%%k==0) factor2=k
			qj=level_vec[factor2]
			freq_levels=matrix(-n/(qi*qj),qi,qj) #freq_levels is the matrix to contain freq of level pairs in the columns i and j
			for (i1 in 1:n)
			{
				freq_levels[D[i1,i],D[i1,j]]=freq_levels[D[i1,i],D[i1,j]]+1
				freq_levels1[D[i1,i],D[i1,i]]=freq_levels1[D[i1,i],D[i1,i]]+1
			}
			for (a in  1:qi)
			{
				for (b in 1:qj)
				{
					freq_levels[a,b]=(freq_levels[a,b])**2										
				}
				for (b in 1:qi)
				{
					freq_levels1[a,b]=(freq_levels1[a,b])**2
				}
			}			
			temp1=sum(freq_levels)	
			temp2=sum(freq_levels1)			
			if (temp1/temp2==1) aliased=aliased + 1 #checks the aliasness of two columns i and j and counts number of such pairs					
		}
	}	
	aliased
}

LB_fNOD=function(m,n,level_vec,k)
{
	#The function returns lower bound of EfNOD criteria for a mixed level design
	sum1=0
	sum2=0
	for (i in 1:(m-1))
	{
		factor1=i%%k
		if (i%%k==0) factor1=k
		qi=level_vec[factor1]
		sum1=sum1+1/qi
		for (j in (i+1):m)
		{
			factor2=j%%k
			if (j%%k==0) factor2=k
			qj=level_vec[factor2]
			sum2=sum2+1/(qi*qj)
		}
	}	
	sum1=sum1+1/level_vec[k]
	sum2=2*sum2
	C=n*m/(m-1)-n*n*(sum1+sum2)/(m*(m-1))
	psi=(n*sum1-m)/(n-1)
	gama=floor(psi)	
	L_fNOD=n*(n-1)*((gama+1-psi)*(psi-gama)+psi**2)/(m*(m-1))+C
	return(L_fNOD)
}

Design_swap=function(D,m,n,k,i,j)
{
#The function generates the design after swaping the column i and j of the generator vector
	t_i=D[1,i]
	t_j=D[1,j]		
	for (row in 1:(n-1))
	{
		coli_plusk=i+(row-1)*k
		colj_plusk=j+(row-1)*k
		if ((coli_plusk)>m) coli_plusk= (coli_plusk) %% m
		if ((colj_plusk)>m) colj_plusk= (colj_plusk) %% m
		D[row,coli_plusk]=t_j
		D[row,colj_plusk]=t_i		
	}
D
}	

mxkssd=function(m,n,level_vec,k,mef)
{
#The function tries to obtain efficient multilevel supersaturated design that has efficiency more than 'mef' by interchanging factors of the generator vector. Here, first position of
#the generator vector is interchanges with most important factor, then second position with most important factor and so on.
#After full iteration, design may not have chisquare efficiency 1. In that case recall the function. mef is the minimum efficiency required, should be between 0 to 1.
	#time of execution
	stime=proc.time()
	L_fNOD=LB_fNOD(m,n,level_vec,k)
	aliased=100000
	trial=0
	while (aliased>0 && trial<=100)  #aliased=1 means two columns in the design are fully aliased
	{
		trial=trial+1
		aliased=0
		gvec=genvec(m,n,level_vec,k)
		D=D_from_G(gvec,m,n,k)	
		outersuccess=1
		final_success=0 #final_success=1 implies a required design was found
		while (outersuccess!=0)  #outersuccess=1 implies that after a full round of iteration the criterion value decreased and a better design was found
		{					
			outersuccess=0
			outer_EfNOD=fNOD(D,m,n,level_vec,k)[[1]]
			loop=0			
			total=k*(n-1)*(n-2)/2
			pb <- txtProgressBar(min = 0, max = total, style=3)	#progress bar variable
			#pb <- winProgressBar(min = 0, max = total)
			i=1				
			while (i<=(m-1) && (final_success==0))
			{
				EfNOD=fNOD(D,m,n,level_vec,k)[[1]]
				min_EfNOD=EfNOD
				if (mef*EfNOD<=L_fNOD) final_success=1
				D1=D								
				j=i+k			
				while ((j<=m) && (final_success==0))
				{
					loop=loop+1 #loop is used for creating progress bar
					Sys.sleep(0.1)
   					#setWinProgressBar(pb, loop) #creates progress bar
					setTxtProgressBar(pb, loop) #creates progress bar
					if (D[1,i]!=D[1,j])
					{
						D_temp=Design_swap(D,m,n,k,i,j)
						EfNOD_temp=fNOD(D_temp,m,n,level_vec,k)[[1]]					
						if (EfNOD_temp<min_EfNOD) 
						{
							min_EfNOD=EfNOD_temp
							D1=D_temp
							success=1							
						}
						if (mef*min_EfNOD<=L_fNOD) final_success=1 #for terminating the outer while loop
					}
					j=j+k				
				}
				D=D1
				i=i+1			
			}
			close(pb) #closes progress bar
			if ((min_EfNOD<outer_EfNOD) & (final_success==0))
			{
				outersuccess=1
				outer_EfNOD=min_EfNOD
			}			 		
		}
		aliased=alias(D,m,n,level_vec,k)	
	}	
	Deff=L_fNOD/min_EfNOD	
	max_fNOD=fNOD(D,m,n,level_vec,k)[[2]]	
	t_taken=proc.time()-stime  #t_taken is time taken to generate the design
	genv=D[1,]
	result=list(m=m,n=n,level_vector=level_vec,k=k,generator.vector=genv,design=D,EfNOD.efficiency=Deff,max.fNOD=max_fNOD,time.taken=t_taken, number.aliased.pairs=aliased)	
	return(result)	
}