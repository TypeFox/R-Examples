emstat.norm <-
function(x,outnull,C,len,niter,tol,k)
#This function computes the EM-test statistics for normal mixture. 
#x:        data; can be either a vector or a matrix with the 1st column being the observed values 
#          		and the 2nd column being the corresponding frequency. 
#outnull:  output from phi0.pois function. 
#C:        optional tuning parameter for EM-test procedure. 
#len:      number of initial values chosen for the EM-algorithm.
#niter:    least number of iterations for all initial values in the EM-algorithm.
#tol:      tolerance value for the convergence of the EM-algorithm. 
#k:	     number of EM iterations to obtain EM-test statistic.
{
	theta0=outnull$theta
	m0=length(theta0)	

	###Calculate eta_h's
	eta=rep(0,m0+1)
	eta[1]=min(x)
	eta[m0+1]=max(x)
	if(m0>1)
	{
		for(i in 2:m0)
			eta[i]=(theta0[i-1]+theta0[i])/2
	}

	###Calculate the collection of beta	
	bbeta=c()
	for(h in 1:m0) 
	{
		bbeta=rbind(cbind(bbeta,rep(0.1,3^{h-1})),
		cbind(bbeta, rep(0.3,3^{h-1})),
		cbind(bbeta, rep(0.5,3^{h-1})))
	}

	pln1=c()
	para=c()
	##For each beta, calculate the statistic m_n^{(k)}
	for(j in 1:(3^m0))
	{
		output=c()
		beta=bbeta[j,]
		para0=maxmm.norm0(x,beta,theta0,len,niter,tol)
		alpha=para0$alpha
		theta1=para0$theta1
		theta2=para0$theta2
		
		###Iteration starts from here###
		for(i in 1:k)
		{
			###E-step of EM-algorithm
			alpha1=alpha*beta	
			alpha2=alpha*(1-beta)
			pdf.part1=apply(as.matrix(theta1,ncol=1),1,dnorm,x=x,sd=1)
			pdf.part2=apply(as.matrix(theta2,ncol=1),1,dnorm,x=x,sd=1)
			pdf.component1=t(t(pdf.part1)*alpha1)+1e-100/m0
			pdf.component2=t(t(pdf.part2)*alpha2)+1e-100/m0
			pdf=apply(pdf.component1,1,sum)+apply(pdf.component2,1,sum)	
			w1=pdf.component1/pdf	
			w2=pdf.component2/pdf
			
			###M-step of EM-algorithm
			for(h in 1:m0)
			{
				if (sum(w1[,h])/(sum(w1[,h])+sum(w2[,h]))<=0.5)
					beta[h]=min((sum(w1[,h])+C)/(sum(w1[,h])+sum(w2[,h])+C),0.5)
				else
					beta[h]=max((sum(w1[,h]))/(sum(w1[,h])+sum(w2[,h])+C),0.5)
			}
			alpha=apply(w1+w2,2,mean)
			alpha1=alpha*beta
			alpha2=alpha*(1-beta)
			theta1=apply(w1*x,2,sum)/apply(w1,2,sum)
			theta2=apply(w2*x,2,sum)/apply(w2,2,sum)
		}		
		###Compute the penalized log-likelihood value and EM-test statistic
		pdf.part1=apply(as.matrix(theta1,ncol=1),1,dnorm,x=x,sd=1)
		pdf.part2=apply(as.matrix(theta2,ncol=1),1,dnorm,x=x,sd=1)
		pdf.component1=t(t(pdf.part1)*alpha1)+1e-100/m0	
		pdf.component2=t(t(pdf.part2)*alpha2)+1e-100/m0
		pdf=apply(pdf.component1,1,sum)+apply(pdf.component2,1,sum)
		para=rbind(para,c(alpha1,alpha2,theta1,theta2))
		pln1[j]=sum(log(pdf))+C*sum(log(1-abs(1-2*beta)))
	}
	index=which.max(pln1)
	para1=para[index,]
	emnk=2*(pln1[index]-outnull$loglik)

	if(emnk<0) 
		emnk=0
	c(emnk,para1)
}
