emtest.binom <-
function(x,size,m0=1,C=NULL,inival=NULL,len=10,niter=50,tol=1e-6,k=3,rformat=FALSE)
#This function compute the EM-test statistic and p-value for the hypothesis H0:m=m0.
#x: 		data, can be either a vector or a matrix with the 1st column being the observed values 
#       		and the 2nd column being the corresponding frequencies. 
#size:  	number of trials.
#m0: 		order under null hypothesis.
#C: 		optional tuning parameter for EM-test procedure.
#inival:	initial values chosen for the EM-algorithm to compute the PMLE under the null model 
#len: 	number of initial values chosen for the EM-algorithm.
#niter:     least number of iterations for all initial values in the EM-algorithm.
#tol: 	tolerance value for the convergence of the EM-algorithm.
#k: 		number of EM iterations to obtain EM-test statistic.
#rformat	format for output, rformat=T means the format of output is determined by R software.
#		rformat=F means the format of output is determined by our default setting. When the output is
#		larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
#		it is determined by signif(output,3).
{
	if (is.data.frame(x))
	{	
		if (ncol(x)==2)
			x=as.matrix(x)
		if (ncol(x)==1 | ncol(x)>2)
			x=x[,1]
	}
	if(is.vector(x))
	{
		y=as.matrix(table(x))
		count=as.numeric(rownames(y))
		freq=y[,1]
		x=cbind(count,freq)
	}
	freq=as.numeric(x[,2])
	n=sum(freq)	
	
	###MLE of parameters under the null model	
	outnull=phi0.binom(x,size,m0,0,inival,len,niter,tol)
	alpha=outnull$alpha
	theta=outnull$theta
	t0=rbind(alpha,theta)

	###Size of penalized function
	if (m0==1)
	{
		if (is.null(C))
			C=0.54
		ah=c(0.5,0.5)
	}
	if (m0==2)
	{
		tb2=tb2.binom(outnull$alpha,outnull$theta,size)
		if (is.null(C))
			C=0.5*exp(5-10.6*tb2[1,2]-123/n)/(1+exp(5-10.6*tb2[1,2]-123/n) )
		ah=c(0.5-acos(tb2[1,2])/2/pi,0.5,acos(tb2[1,2])/2/pi)
	}
	if (m0==3)
	{
		tb2=tb2.binom(outnull$alpha,outnull$theta,size)
		if (is.null(C))
			C=0.5*exp(3.3-5.5*tb2[1,2]-5.5*tb2[2,3]-165/n)/(1+exp(3.3-5.5*tb2[1,2]-5.5*tb2[2,3]-165/n)) 
		
		a0=0.5-acos(tb2[1,2])/4/pi-acos(tb2[1,3])/4/pi-acos(tb2[2,3])/4/pi
		a2=0.5-a0
		w123=(tb2[1,2]-tb2[1,3]*tb2[2,3])/sqrt(1-tb2[1,3]^2)/sqrt(1-tb2[2,3]^2)
		w132=(tb2[1,3]-tb2[1,2]*tb2[3,2])/sqrt(1-tb2[1,2]^2)/sqrt(1-tb2[3,2]^2)
		w231=(tb2[2,3]-tb2[2,1]*tb2[3,1])/sqrt(1-tb2[2,1]^2)/sqrt(1-tb2[3,1]^2)
		a1=0.75-acos(w123)/4/pi-acos(w132)/4/pi-acos(w231)/4/pi
		a3=0.5-a1
		ah=c(a0,a1,a2,a3)
	}
	if (m0>3)
	{
		if (is.null(C))
			C=0.5
		tb2=tb2.binom(outnull$alpha,outnull$theta,size)
		ah=emtest.thm3(tb2,N=10000,tol=1e-8)
	}

	out=emstat.binom(x,outnull,size,C,len,niter,tol,k)
	emnk=out[1]
	alpha=out[2:(2*m0+1)]
	theta=out[(2*m0+2):(4*m0+1)]
	t1=rbind(alpha,theta)
	p=sum(ah*pchisq(emnk,0:m0,lower.tail = F))

	if (rformat==F)
	{
		t0=rousignif(t0)
		t1=rousignif(t1)
		emnk=rousignif(emnk)
		p=rousignif(p)
		C=rousignif(C)
	}

	list('MLE of parameters under null hypothesis (order = m0)'=t0,
	'Parameter estimates under the order = 2m0'=t1,
	'EM-test Statistic'=emnk,
	'P-value'=p,
	'Level of penalty'=C)
}
