phi0.exp <-
function(x,m0,lambda,inival,len,niter,tol)
#This function computes the PMLE of mixing distribution under null hypothesis for exponential mixture. 
#x:      data; can be either a vector or a matrix with the 1st column being the observed values 
#        		and the 2nd column being the corresponding frequencies. 
#m0:     order under null hypothesis.
#lambda:	size of penalty function of mixing proportions.
#inival: initial values for the EM-algorithm.
#len: 	number of initial values chosen for the EM-algorithm.
#niter:     least number of iterations for all initial values in the EM-algorithm.
#tol: 	tolerance value for the convergence of the EM-algorithm.
{
	n=length(x)
	if(m0>1)
	{
		if (is.vector(inival))
			inival=t(inival)
		if (is.null(inival)==F)
			len=nrow(inival)
		output=c()
		for(i in 1:len)
		{
			if (is.null(inival))
			{
				alpha=runif(m0,0,1)
				alpha=alpha/sum(alpha)
				theta=sort(runif(m0,min(x),max(x)))
			}
			else
			{
				alpha=inival[i,1:m0]
				alpha=alpha/sum(alpha)
				theta=sort(inival[i,(m0+1):(2*m0)])
			}
			for (j in 1:niter)###run niter EM-iterations first
			{
				pdf.component=t(t(apply(as.matrix(1/theta,ncol=1),1,dexp,x=x))*alpha)+1e-100/m0
				pdf=apply(pdf.component,1,sum)
				w=pdf.component/pdf
				alpha=(apply(w,2,sum)+lambda)/(n+m0*lambda)
				theta=apply(w*x,2,sum)/apply(w,2,sum)
			}
			pdf.component=t(t(apply(as.matrix(1/theta,ncol=1),1,dexp,x=x))*alpha)+1e-100/m0
			pdf=apply(pdf.component,1,sum)
			pln=sum(log(pdf))+lambda*sum(log(alpha))
			output=rbind(output,c(alpha,theta,pln))
		}
		index=which.max(output[,(2*m0+1)])
		alpha=output[index,1:m0]
		theta=output[index,(m0+1):(2*m0)]
		pln0=output[index,(2*m0+1)]
		err=1
		t=0
		pdf.component=t(t(apply(as.matrix(1/theta,ncol=1),1,dexp,x=x))*alpha)+1e-100/m0
		pdf=apply(pdf.component,1,sum)
		while(err>tol & t<2000)###EM-iteration with the initial value with the largest penalized log-likelihood
		{
			w=pdf.component/pdf
			alpha=(apply(w,2,sum)+lambda)/(n+m0*lambda)
			theta=apply(w*x,2,sum)/apply(w,2,sum)
			pdf.component=t(t(apply(as.matrix(1/theta,ncol=1),1,dexp,x=x))*alpha)+1e-100/m0
			pdf=apply(pdf.component,1,sum)
			pln1=sum(log(pdf))+lambda*sum(log(alpha))
			err=abs(pln1-pln0)
			pln0=pln1
			t=t+1
		}
		ln=pln1-lambda*sum(log(alpha))
		index=sort(theta,index.return=TRUE)$ix
		alpha0=alpha[index]
		theta0=theta[index]

		list("alpha"=alpha0,"theta"=theta0,"loglik"=ln,"ploglik"=pln1)
	}
	else 
	{
		theta0=mean(x)
		list("alpha"=1,"theta"=theta0,"loglik"=sum(log(dexp(x,rate=1/theta0))),
			"ploglik"=sum(log(dexp(x,rate=1/theta0))))
	}
}
