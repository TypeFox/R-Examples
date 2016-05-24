is.wholenumber=function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
library(lpSolve)
ccmat_LP=function(v,b,k)
{
	if (v<=k) stop ("v should be greater than k") 
	r=floor(b*k/v)
	l1=v-(b*k-v*r)
	l2=v-l1
	lambda=floor(r*(k-1)/(v-1))
	n12=r*(k-1)-lambda*(v-1)
	n11=v-1-n12
	if (n11>=(k-1))
	{
		n21=n11-k+1
		n22=v+k-n11-2
	} 
	if (n11<(k-1))
	{
		n21=v-k+n11
		n22=k-1-n11
	}
	
	if ((b*k/v==r) & (r*(k-1)/(v-1)==lambda))
	{
		NNP=matrix(lambda,v,v) 
		for (i in 1:l1) NNP[i,i]=r
		stop(return(NNP))
	}
	obj=matrix(0,v,v)	
	dim(obj)=c(1,v*v)
	constr1=NULL
	for (i in 1:v) 
	{
		temp=matrix(0,v,v)
		temp[i,]=1
		temp[i,i]=1-k
		temp=t(temp)
		dim(temp)=c(1,v*v)
		constr1=rbind(constr1,temp)
	}
	constr2=NULL
	for (i in 1:(v-1))
	{
		for (j in (i+1):v)
		{
			temp=matrix(0,v,v)
			temp[i,j]=1
			temp[j,i]=-1
			temp=t(temp)
			dim(temp)=c(1,v*v)
			constr2=rbind(constr2,temp)
		}	
	}
	constr3=matrix(1,v,v)
	for (i in 1:v) constr3[i,i]=0
	constr3=t(constr3)
	dim(constr3)=c(1,v*v)
	constr4=NULL
	for (i in 1:v) 
	{
		temp=matrix(0,v,v)
		temp[i,i]=1				
		dim(temp)=c(1,v*v)
		constr4=rbind(constr4,temp)
	}	
	constr5=NULL
	for (i in 1:(v-1))
	{
		for (j in (i+1):v)
		{
			temp=matrix(0,v,v)
			temp[i,j]=1
			temp=t(temp)
			dim(temp)=c(1,v*v)
			constr5=rbind(constr5,temp)
		}	
	}
	constr51=constr5
	constr6=NULL
	for (i in 1:(v-1))
	{
		for (j in (i+1):v)
		{
			temp=matrix(0,v,v)
			temp[i,i]=1
			temp[i,j]=-1
			temp=t(temp)
			dim(temp)=c(1,v*v)
			constr6=rbind(constr6,temp)			
		}					
	}	
	constr=rbind(constr1,constr2,constr3,constr4,constr5,constr51)
	dir1=rep("=", times=(v))
	dim(dir1)=c(v,1)	
	dir2=rep("=",times=(v*(v-1)/2))	
	dim(dir2)=c(v*(v-1)/2,1)	
	dir3=rep("=",times=(1)) 
	dir4=rep("==", times=(v))
	dim(dir4)=c(v,1)
	dir5=rep(">=",times=(v*(v-1)/2))	
	dim(dir5)=c(v*(v-1)/2,1)
	dir51=rep("<=",times=(v*(v-1)/2))	
	dim(dir51)=c(v*(v-1)/2,1)
	dir6=rep(">",times=(v*(v-1)/2))	
	dim(dir6)=c(v*(v-1)/2,1)						
	dir=rbind(dir1,dir2,dir3,dir4,dir5,dir51)	
	rhs1=matrix(0,v,1)
	rhs2=matrix(0,v*(v-1)/2,1)
	rhs3=b*k*(k-1)	
	rhs4=matrix(r,l1,1)
	if (l2>0) rhs4=rbind(rhs4,matrix((r+1),l2,1))
	rhs5=matrix(floor(r*(k-1)/(v-1)),v*(v-1)/2,1)
	if (l2==0 | (l2>0 & n11>=(k-1))) rhs51=matrix((floor(r*(k-1)/(v-1))+1),v*(v-1)/2,1) 
	if (l2>0 & n11<(k-1)) rhs51=rbind(matrix((floor(r*(k-1)/(v-1))+1),(v*(v-1)/2-l2*(l2-1)/2),1),matrix((floor(r*(k-1)/(v-1))+2),l2*(l2-1)/2,1))
	rhs6=matrix(0,v*(v-1)/2,1)
	rhs=rbind(rhs1,rhs2,rhs3,rhs4,rhs5,rhs51)	
	sol=lp (direction = "min", obj, constr, dir, rhs,transpose.constraints = TRUE, all.int=TRUE)	
	if (sol[[28]]==0) 
	{		
		row=sol[[12]]
		dim(row)=c(v,v)		
	} else {
			rhs51=matrix((floor(r*(k-1)/(v-1))+3),(v*(v-1)/2),1)
			rhs=rbind(rhs1,rhs2,rhs3,rhs4,rhs5,rhs51)	
			sol=lp (direction = "min", obj, constr, dir, rhs,transpose.constraints = TRUE, all.int=TRUE)	
			if (sol[[28]]==0) 
			{		
				row=sol[[12]]
				dim(row)=c(v,v)		
			} else return(0)
		}	
	return(row)	
}
#######################################################################
ccmat=function(v,b,k)
{
	x=NULL	
	lambda=floor(b*k*(k-1)/(v*(v-1)))
	n2=b*k*(k-1)/2-v*(v-1)*lambda/2
	n1=v*(v-1)/2-n2	
	r=floor(b*k/v)
	l2=b*k-v*r 
	l1=v-l2
	n12=r*(k-1)-lambda*(v-1)
	n11=v-1-n12
	n22=(r+1)*(k-1)-lambda*(v-1)
	n21=v-1-n22
	if ((((b*k/v)-r)==0) & (((b*k*(k-1)/(v*(v-1)))-lambda)==0)) 
	{
		NNP=matrix(lambda,v,v) 
		for (i in 1:l1) NNP[i,i]=r
		stop(return(NNP))
	} else NNP=matrix((lambda+1),v,v)	
	for (i in 1:l1) NNP[i,i]=r
	if (l2>0) 
	{
		for (i in (l1+1):v) NNP[i,i]=r+1
	}
	parm=c(v,b,k,r,lambda,l1,l2,n11,n12,n21,n22)
	if (all(parm>=0))
	{	
		row_cnt=matrix(0,l1,1)
		col_cnt=matrix(0,l1,1)
		rem_row=c(1:l1)
		rem_col=c(1:l1)
		trial1=0
		while (length(rem_row)!=0 & trial1<=10000)
		{
			row=sample(rem_row,1)
			col=sample(rem_col,1)
			if (row!=col & NNP[row,col]==(lambda+1) & row_cnt[row,]<n11 & row_cnt[col,]<n11)
			{
				NNP[row,col]=lambda
				NNP[col,row]=lambda
				row_cnt[row,]=row_cnt[row,]+1
				col_cnt[col,]=col_cnt[col,]+1
				row_cnt[col,]=row_cnt[col,]+1
				col_cnt[row,]=col_cnt[row,]+1
				if (row_cnt[row,]==n11) rem_row=setdiff(rem_row,row)
				if (row_cnt[col,]==n11) rem_row=setdiff(rem_row,col)
				if (col_cnt[col,]==n11) rem_col=setdiff(rem_col,col)
				if (col_cnt[row,]==n11) rem_col=setdiff(rem_col,row)
			}
			trial1=trial1+1
		}
		if (l2>0)
		{
			row_cnt=matrix(0,l2,1)
			col_cnt=matrix(0,l2,1)
			rem_row=c((l1+1):v)
			rem_col=c((l1+1):v)
			trial2=0
			while (length(rem_row)!=0 && trial2<=10000)
			{
				length_rem_row=length(rem_row)
				length_rem_col=length(rem_col)
				row=rem_row[sample(length_rem_row,1)]
				col=rem_col[sample(length_rem_col,1)]
				if (row!=col && NNP[row,col]==(lambda+1) && row_cnt[row-l1,]<n21 && row_cnt[col-l1,]<n21)
				{
					NNP[row,col]=lambda
					NNP[col,row]=lambda
					row_cnt[row-l1,]=row_cnt[row-l1,]+1
					col_cnt[col-l1,]=col_cnt[col-l1,]+1
					row_cnt[col-l1,]=row_cnt[col-l1,]+1
					col_cnt[row-l1,]=col_cnt[row-l1,]+1
					if (row_cnt[row-l1,]==n21) rem_row=setdiff(rem_row,row)
					if (row_cnt[col-l1,]==n21) rem_row=setdiff(rem_row,col)
					if (col_cnt[col-l1,]==n21) rem_col=setdiff(rem_col,col)
					if (col_cnt[row-l1,]==n21) rem_col=setdiff(rem_col,row)
				}
				trial2=trial2+1
			}
		}
		if(trial1<10000) return(NNP) else return(matrix(0,v,v))
	} else return(parm)	
}
##################################################################################################
N_to_design=function(N)
{
	design=NULL
	v=nrow(N)
	b=ncol(N)
	kvec=t(N)%*%matrix(1,v,1)
	k=max(kvec)	
	for (i in 1:b)
	{
		temp=which(N[,i]>0)
		if (length(temp)<k) padzero=rep(0,times=(k-length(temp))) else padzero=NULL
		design=rbind(design,c(temp,padzero))
	}
	return(design)
}
##################################################################################################
is.proper=function(N)
{
	v=nrow(N)
	kvec=t(N)%*%matrix(1,v,1)
	k=sum(N[,1])
	if (all(kvec==k)) proper=1 else proper=0
	return(proper)
}
##################################################################################################
Cmatrix=function(N)
{
	v=nrow(N)	
	NNP=N%*%t(N)	
	R=matrix(0,v,v);
	for (i in 1:v)
	{
		R[i,i]=NNP[i,i];
	}		
	if (is.proper(N)==1) 
	{
		k=sum(N[,1])
		C=R-NNP/k
	} else {
			kvec=t(N)%*%matrix(1,v,1)
			kvec=c(kvec)
			b=length(kvec)
			Kinv=diag(1/kvec,nrow=b,ncol=b)
			C=R-N%*%Kinv%*%t(N)	
		}				
	return(C)
}
##################################################################################################
is.connected=function(N)
{
	v=nrow(N)
	C=Cmatrix(N)
	dt=det(C+matrix(1/v,v,v))
	dt=round(dt,digits=6)
	if (dt!=0)  connected=1 else connected=0 
	return(connected)
}
##################################################################################################
A_eff.NNP=function(b,k,NNP)
{
	v=nrow(NNP)		
	R=diag(diag(NNP))		
	C=R-NNP/k	
	dt=det(C+matrix(1/v,v,v))
	dt=round(dt,digits=6)
	if (dt!=0) 
	{
		Eigen_out=eigen(C)
		eig_values=Eigen_out$values
		eig_values=round(eig_values, digits=6)
		sum=0
		for (i in 1:v)
		{
			if ((eig_values[i]) !=0)  sum=sum+(1/(eig_values[i]))
		}	
		LB_Aeff=((v-1)**2)/(b*(k-1))/sum		
	} else LB_Aeff=0
	return(LB_Aeff)
	
}
##################################################################################################
A_eff=function(N)
{
	if (is.proper(N)==1) k=sum(N[,1]) else stop ("design is not  proper")
	v=nrow(N)
	b=ncol(N)	
	C=Cmatrix(N)
	Eigen_out=eigen(C)
	eig_values=Eigen_out$values
	eig_values=round(eig_values, digits=6)  # to consider decimal places upto 2
	sum=0
	for (i in 1:v)
	{
		if ((eig_values[i]) !=0)  sum=sum+(1/(eig_values[i]))  #in place of zero eigen values, we have used 0.00000001
	}	
	LB_Aeff=((v-1)**2)/(b*(k-1))/sum
	return(LB_Aeff)	
}
##################################################################################################
D_eff=function(N)
{
	if (is.proper(N)==1) k=sum(N[,1]) else stop ("design is not  proper")
	v=nrow(N)
	b=ncol(N)	
	C=Cmatrix(N)
	Eigen_out=eigen(C)
	eig_values=Eigen_out$values
	eig_values=round(eig_values, digits=6)  # to consider decimal places upto 2
	product=1
	for (i in 1:v)
	{
		if ((eig_values[i])!=0)  product=product*(1/(eig_values[i]))
	}
	GM=product**(1/(v-1))	
	LB_Deff=(v-1)/(b*(k-1)*GM)	
	return(LB_Deff)
}
##################################################################################################
check.validity.NNP=function(NNP,k)
{
	v=nrow(NNP)
	temp=0
	for (i in 1:v)
	{
		if (is.wholenumber(NNP[i,i]))
		{
			if (((NNP[i,i])*(k-1))==(sum(NNP[i,])-NNP[i,i])) temp=temp+1
		}
	}
	if (temp==v) valid=1 else valid=0
	return(valid)
}
#################################################################################################
interchange.NNP=function(b,k,NNP)
{
	Aeff=A_eff.NNP(b,k,NNP)
	v=nrow(NNP)
	i1=1	
	while (i1<=(v-1))
	{
		j1=i1+1
		flag=0
		while (j1<=v & flag==0)
		{
			i2=i1
			while (i2<=(v-1) & flag==0)
			{
				j2=i2+1
				while (j2<=v & flag==0)
				{	
					if ((i1!=i2) | (j1!=j2)) 
					{
						if (NNP[i1,j1]!=NNP[i2,j2])
						{
							NNP.temp=NNP
							NNP.temp[i1,j1]=NNP[i2,j2]
							NNP.temp[i2,j2]=NNP[i1,j1]
							NNP.temp[j1,i1]=NNP.temp[i1,j1]
							NNP.temp[j2,i2]=NNP.temp[i2,j2]
							NNP.temp[i1,i1]=(sum(NNP.temp[i1,])-NNP[i1,i1])/(k-1)
							NNP.temp[j1,j1]=(sum(NNP.temp[j1,])-NNP[j1,j1])/(k-1)
							NNP.temp[i2,i2]=(sum(NNP.temp[i2,])-NNP[i2,i2])/(k-1)
							NNP.temp[j2,j2]=(sum(NNP.temp[j2,])-NNP[j2,j2])/(k-1)
							if (check.validity.NNP(NNP.temp,k)==1) 
							{
								Aeff.temp=A_eff.NNP(b,k,NNP.temp)
								if (Aeff.temp>Aeff) 
								{
									Aeff=Aeff.temp
									NNP=NNP.temp
									i1=0
									flag=1
								}
							}
						}
					}
					j2=j2+1	
				}
				i2=i2+1
			}
			j1=j1+1
		}
		i1=i1+1
		
	}
	return(NNP)
}
################################################################
LIP=function(v,b,kvec,NNPo,N1,T,rownum,relaxed)
{
	kvec_obt=t(N1)%*%matrix(1,nrow(N1),1)
	w=matrix(0,1,b)
	for (j in 1:b)
	{
		if (kvec_obt[j,]==0) w[,j]=1
		else w[,j]=1/kvec_obt[j,]
	}		
	obj=w	#kvec_obt is the kvec obtained till (i-1) step
	constr1=matrix(1,1,b)		
	constr2=matrix(0,b,b)
	for (j in 1:b)
	{
		constr2[j,j]=1	
	}
	constr3=N1 	
	constr4=T
	constr=rbind(constr1,constr2,constr3,constr4)
	if (relaxed>0)  
	{
		constr=rbind(constr1,constr2,constr4)
	}
	dir1=rep("=", times=(1))	
	dir2=rep("<=",times=(b))	
	dim(dir2)=c(b,1)	
	dir3=rep("=",times=(nrow(N1)))  
	dim(dir3)=c(nrow(N1),1)	
	dir4=rep("<",times=(nrow(constr4)))
	dim(dir4)=c(nrow(constr4),1)					
	dir=rbind(dir1,dir2,dir3,dir4)	
	if (relaxed>0) 
	{
		dir=rbind(dir1,dir2,dir4)
	}
	rhs1=NNPo[rownum,rownum]
	rhs2=kvec-kvec_obt 
	rhs3=matrix(0,nrow(N1),1)
	for (j in 1:nrow(N1))
	{
		if (sum(N1[j,])>0) rhs3[j,]=NNPo[j,rownum] else rhs3[j,]=0
	}		
	rhs4=matrix((NNPo[rownum,rownum]-0.5),nrow(constr4),1)
	rhs=rbind(rhs1,rhs2,rhs3,rhs4)
	if (relaxed>0) 
	{
		rhs=rbind(rhs1,rhs2,rhs4)
	}
	types=rep("B", times=b)	
	sol=lp (direction = "max", obj, constr, dir, rhs,transpose.constraints = TRUE, all.bin=TRUE, use.rw=TRUE)
	if (sol[[28]]==0) 
	{		
		row=sol[[12]]
		dim(row)=c(1,b)
		if (rownum>nrow(N1)) N1=rbind(N1,row) else N1[rownum,]=row		
	} 	
	return(N1)	
}
##################################################################################################
detect=function(v,b,kvec,NNPo,N1,T,relaxed)
{
	row_detected=0	
	result=0
	k0=1   
	#while (k0<=nrow(N1) & row_detected==0)
	while (k0<=min(4,nrow(N1)) & row_detected==0)
	{
		row_indices=combn(nrow(N1),k0)
		nr=ncol(row_indices)
		j=1
		while(j<=nr & row_detected==0)
		{
			rows=row_indices[,j]
			T_temp=rbind(T,N1[rows,])			
			N1_temp=N1
			N1_temp[rows,]=matrix(0,1,b)
			cnt=0
			for (m in 1:k0)
			{
				rownum=rows[m]				
				N1_temp=LIP(v,b,kvec,NNPo,N1_temp,T_temp,rownum,relaxed)
				if (sum(N1_temp[rownum,])>0) cnt=cnt+1
			}
			if (cnt==k0) {
					row_detected=1
					result=list(rows,N1_temp)
				     }
			j=j+1	
		}
		k0=k0+1				
	}
	return(result)
}
##################################################################################################
ibdgen=function(v,b,k,NNPo,ntrial,pbar)
{
	R=matrix(0,v,v);
	for (i in 1:v)
	{
		R[i,i]=NNPo[i,i];
	}		
	C=R-NNPo/k	
	dt=det(C+matrix(1/v,v,v))
	dt=round(dt,digits=6)
	if (dt!=0)  connected=1 else connected=0
	if (connected==1)
	{	
		kvec=matrix(k,b,1)
		trial=1
		success=0				
		while(trial<=ntrial & success==0)
		{
			if (pbar==TRUE) 
			{
				if (Sys.info()[[1]]=="Windows") pb = winProgressBar(title = "progress bar", min = 0, max = v, width = 400) else pb=txtProgressBar(min = 0, max = v, style=3)
			}
			N1=matrix(0,1,b)
			col=sample(b,(NNPo[1,1]))
			N1[1,col]=1
			T=matrix(0,1,b)
			i=2	
			decision=0
			relaxed=0				
			while (i<=v & decision==0)
			{
				nt=nrow(T)
				if (nt>5*v)  
				{
					T=matrix(0,1,b)
					decision=1 						
				}
				N1=LIP(v,b,kvec,NNPo,N1,T,i,relaxed)
				if (nrow(N1)<i ) 
				{
					temp=detect(v,b,kvec,NNPo,N1,T,relaxed)
					rows=temp[[1]]
					if (all(rows>0)) 
					{
						T=rbind(T,N1[rows,])
						N1=temp[[2]]												
					} else {
							decision=1 								
					           }			
				} 					
				if (nrow(N1)<i & trial==ntrial)
				{
					relaxed=1							
				} 
				Sys.sleep(0.1)
				if (pbar==TRUE) 
				{
					if (Sys.info()[[1]]=="Windows") setWinProgressBar(pb, i,title=paste(round((i-1)*100/v, 0),"% done,","row=",i, ",trial=",trial, ",tabulist=", nt)) else  setTxtProgressBar(pb, i)
				}
				i=nrow(N1)+1	
			}
			if (pbar==TRUE) close(pb)
			trial=trial+1			
			if (nrow(N1)==v) 
			{
				connected1=is.connected(N1)
				if (connected1==1)
				{
					success=1
					result=N1
				} else {
					design="Connected design not found"	
					result=list(v=v,b=b,k=k,design=design)
	         		       }
			} else {
				design="Design not found"	
				result=list(v=v,b=b,k=k,design=design)
	         	       }
		}
	} else {
		design="Suitable cocurrence matrix of a connected design was not found"	
		result=list(v=v,b=b,k=k,design=design)
	       }	
	return(result)			
}
##################################################################################################
do.exist.NBIB=function(v,b,k)
{
	r=floor(b*k/v)
	l1=v-(b*k-v*r)	
	lambda=floor(r*(k-1)/(v-1))
	n11=v-1-(r*(k-1)-lambda*(v-1))	
	if (((n11>=(k-1)) & (l1*(n11-l1+1)<=(v-l1)*(n11-k+1))) | ((n11<(k-1)) & ((n11+1)<=l1) & (l1<=(v-k+n11)) & ((n11*l1)%%2==0) & (((v-l1)*(v-k+n11-l1))%%2==0))) exist=1 else exist=0
	return(exist)
}
##################################################################################################
ibd=function(v,b,k,NNPo,ntrial=5,pbar=FALSE)
{
	stime=proc.time()
	if (v<0 | k<0 | b<0) stop("v,b,k should be positive.")
	if (!is.wholenumber(v) | !is.wholenumber(b) | !is.wholenumber(k)) stop("v,b,k should be integers.")
	if (b*k<v+b-1) stop ("Design parameters do not satisfy conditions for even minimal connectedness.")
	if (!is.wholenumber(ntrial)) stop("ntrial should be integer.")
	valid=0
	if (missing(NNPo))
	{
		try=0
		while (valid==0 & try<11)
		{
			if (try<10) NNPo=ccmat_LP(v,b,k) else NNPo=ccmat(v,b,k)
			if (length(NNPo)==v*v & sum(NNPo)>0) valid=check.validity.NNP(NNPo,k)
			try=try+1
			if (valid==1) NNPo=interchange.NNP(b,k,NNPo)
		}
	} else {
		valid=check.validity.NNP(NNPo,k)
		if (valid==0) stop("Check your concurrence matrix.")
	       }	
	if (valid==1) 
	{
		N=ibdgen(v,b,k,NNPo,ntrial,pbar)
		if (is.matrix(N))
		{
			design=N_to_design(N)
			conc_mat=N%*%t(N)
			Aeff=A_eff(N)
			Deff=D_eff(N)
			t.taken=proc.time()-stime
			result=list(v=v,b=b,k=k, NNP=NNPo, N=N, design=design,conc.mat=conc_mat,A.Efficiency=Aeff, D.Efficiency=Deff,time.taken=t.taken) 
		} else result=N
		
	}else {
		design="Suitable cocurrence matrix of a connected design was not found"	
		 result=list(v=v,b=b,k=k,design=design)
	          }
	return(result)
}
##################################################################################################
bibd=function(v,b,r,k,lambda,ntrial=5,pbar=FALSE)
{
	if (v*r==b*k & lambda*(v-1)==r*(k-1))
	{
		NNPo=matrix(lambda,v,v)
		diag(NNPo)=r
		N=ibdgen(v,b,k,NNPo,ntrial,pbar)
		if (is.matrix(N))
		{	
			NNP=N%*%t(N)
			design=N_to_design(N)
			Aeff=A_eff(N)
			Deff=D_eff(N)
			result=list(v=v,b=b,r=r,k=k,lambda=lambda,design=design,N=N,NNP=NNP,Aeff=Aeff,Deff=Deff)
		} else result="design not found"
	} else result="parameters do not satisfy necessary conditions"
	return(result)
}
##################################################################################################
design_to_N=function(design)
{
	v=max(design)
	b=nrow(design)
	k=ncol(design)
	N=matrix(0,v,b)
	for (i in 1:b)
	{
		for (j in 1:k)
		{
			N[design[i,j],i]=N[design[i,j],i]+1
		}
	}
	return(N)
}
##################################################################################################
library(MASS)
A_eff_tc=function(N,v1,v2,b,k)
{
	kby2=floor(k/2)
	min=99999
	for (x in 0:(kby2-1))
	{
		for (z in 0:b)
		{
			if (x==0 & z==0) g=99999 else 
						{
						  	a=v2*(v1-1)**2
							d=v1*(v2-1)
							C=b*x+z
							A= (k*C-v2*(b*x*x+2*x*z+z))/(v1*k)
							B=(b*k*v1*(k-1)-v2*C*(v1*(k-1)+k)+v2*v2*(b*x*x+2*x*z+z))/(v1*k)
							g=1/A+a/B+d/C
						}
			if (g<min) min=g
		}
	}
	LB=min
	rvec=N%*%matrix(1,b,1)		
	v=v1+v2
	M=Cmatrix(N)
	Minv=ginv(M)
	onev2=matrix(1,v2,1)
	onev1=matrix(1,v1,1)
	iv2=matrix(0,v2,v2)
	diag(iv2)=1
	iv1=matrix(0,v1,v1)
	diag(iv1)=1
	P=cbind(onev2%x%iv1,-(iv2%x%onev1))
	T=P%*%Minv%*%t(P)
	temp=0
	for (i in 1:(nrow(T))) temp=temp+T[i,i]
	e=min/temp	
	return(e)			
}
#############################################################################
NLIP=function(v,b,k,NNPo)
{
	B=combn(v,k) 
	order=ncol(B) 	
	obj=matrix(1,1,order)
	pairs=combn(v,2)
	constr1=matrix(0,v*(v-1)/2,order)		
	for (i in 1:(v*(v-1)/2))
	{
		for (j in 1:order)
		{
			if (all(is.element(pairs[,i],B[,j]))) constr1[i,j]=1
		}
	}
	constr2=matrix(1,1,order)
	constr=rbind(constr1,constr2)
	dir1=rep("==", times=(v*(v-1)/2))
	dim(dir1)=c(v*(v-1)/2,1)
	dir2="=="		
	dir=rbind(dir1,dir2)
	rhs1=matrix(0,v*(v-1)/2,1)		
	count=0
	for (i in 1:(v-1))
	{
		for (j in (i+1):v)
		{
			count=count+1			
			rhs1[count,]=NNPo[i,j]		
		}			
	}
	rhs2=b	
	rhs=rbind(rhs1,rhs2)	
	types=rep("B", times= ncol(constr))		
	sol=lp (direction = "min", obj, constr, dir, rhs,transpose.constraints = TRUE, all.bin=TRUE, use.rw=TRUE)	
	if (sol[[28]]==0) 
	{		
		cols=which(sol[[12]]>0)
		result=t(B[,cols])		
	} else result="No feasible solution"		
	return(result)
}
###############################################################################
btib1=function(v,b,r,r0,k,lambda,lambda0)
{
	NNPo=matrix(lambda,v+1,v+1)
	diag(NNPo)=r
	NNPo[v+1,v+1]=r0
	for (i in 1:v)
	{
		NNPo[i,v+1]=lambda0
		NNPo[v+1,i]=lambda0
	}	
	design=NLIP(v+1,b,k,NNPo)
	if (is.matrix(design)) 
	{
		N=design_to_N(design)
		NNP=N%*%t(N)	
		Aeff=A_eff_tc(N,(v-1),1,b,k)
		out=list(v=v,b=b,r=r,r0=r0,k=k,lambda=lambda,lambda0=lambda0,design=design,N=N,NNP=NNP,Aeff=Aeff)
	} else out=design
	return(out)	
}
###############################################################################
btib=function(v,b,r,r0,k,lambda,lambda0,ntrial=5,pbar=FALSE)
{
	NNPo=matrix(lambda,v+1,v+1)
	diag(NNPo)=r
	NNPo[v+1,v+1]=r0
	for (i in 1:v)
	{
		NNPo[i,v+1]=lambda0
		NNPo[v+1,i]=lambda0
	}
	N=ibdgen((v+1),b,k,NNPo,ntrial,pbar)
	if (is.matrix(N))
	{	
		NNP=N%*%t(N)
		design=N_to_design(N)
		Aeff=A_eff_tc(N,v,1,b,k)
		result=list(v=v,b=b,r=r,r0=r0,k=k,lambda=lambda,lambda0=lambda0,design=design,N=N, NNP=NNP,Aeff=Aeff)
	} else result="design not found"	
	return(result)
}
#################################################################################
ibdtvc=function(v1,v2,b,k,NNPo,ntrial=5,pbar=FALSE)
{
	if (v1<0 |v2<0| k<0 | b<0) stop("v1,v2,b,k should be positive")
	if (missing (ntrial)) ntrial=5
	v=v1+v2
	N=ibdgen(v,b,k,NNPo,ntrial,pbar)
	if (is.matrix(N))
	{	
		NNP=N%*%t(N)
		design=N_to_design(N)
		Aeff=A_eff_tc(N,v1,v2,b,k)
		result=list(v1=v1,v2=v2,b=b,k=k,design=design,N=N, NNP=NNP,Aeff=Aeff)
	} else result="design not found"	
	return(result)	
}
##################################################################################################
is.vb=function(N)
{
	v=nrow(N)
	if(is.connected(N)==1)
	{
		C=Cmatrix(N)
		diag.ele=diag(C)
		offdiag.ele=NULL
		for (i in 1:(v-1))
		{
			for (j in (i+1):v)
			{
				offdiag.ele=rbind(offdiag.ele,C[i,j])
			}
		}
		if (all(diag.ele==diag.ele[1]) & all(offdiag.ele==offdiag.ele[1,])) vb=1 else vb=0
	} else vb=0
	return(vb)
}
##################################################################################################
is.equir=function(N)
{
	b=ncol(N)
	rvec=N%*%matrix(1,b,1)
	r=sum(N[1,])
	if (all(rvec==r)) equir=1 else equir=0
	return(equir)
}
##################################################################################################
is.orthogonal=function(N)
{
	v=nrow(N)
	b=ncol(N)
	n=v*b
	rvec=N%*%matrix(1,b,1)
	kvec=t(N)%*%matrix(1,v,1)
	if (all(N==(rvec%*%t(kvec)/n))) orthogonal=1 else orthogonal=0
	return(orthogonal)
}
###################################################################################################
check.orthogonality=function(M)
{
	orthog=1
	for (i in 1:(nrow(M)-1))
	{
		for (j in (i+1):nrow(M))
		{
			if (sum(M[i,]*t(M[j,]))!=0) orthog=0
		}
	}
	return(orthog)
			
}
##################################################################################################
library(lsmeans)
library(car)
aov.ibd=function(formula,specs,data,contrast,joint=FALSE,details=FALSE,sort=TRUE,by=NULL,alpha = 0.05,Letters="ABCDEFGHIJ",...)
{
	lm.obj=lm(formula,data)
	ANOVA.table=Anova(lm.obj,type="III")
	if (details) res=list(lm.obj=lm.obj,ANOVA.table=ANOVA.table) else res=list(ANOVA.table=ANOVA.table)
	if(!missing(specs))
	{
		lsmeans.obj=lsmeans(lm.obj,specs,data=data)
		cld.obj=cld(lsmeans.obj,details,sort,by,alpha,Letters)
		if (details) res=list(lm.obj=lm.obj,ANOVA.table=ANOVA.table,LSMEANS=cld.obj) else res=list(ANOVA.table=ANOVA.table,LSMEANS=cld.obj)
		if(!missing(contrast))
		{
			nr=nrow(contrast)
			contrast.list=vector("list",nr)
			for (i in 1:nr)
			{
				contrast.list[[i]]=contrast[i,]
				names(contrast.list)[[i]]=rownames(contrast)[i]
			}
			contrast.analysis=contrast(lsmeans.obj,contrast.list)
			if (joint) 
			{
				orthogonal=check.orthogonality(contrast)
				if(orthogonal==1) contrast.analysis=lht(lm.obj,contrast.analysis@linfct) else contrast.analysis="contrasts are not orthogonal. Joint test not possible"
			} 
			if (details) res=list(lm.obj=lm.obj,ANOVA.table=ANOVA.table,LSMEANS=cld.obj,contrast.analysis=contrast.analysis) else res=list(ANOVA.table=ANOVA.table,LSMEANS=cld.obj,contrast.analysis=contrast.analysis)
		}
	}
	return(res)
}