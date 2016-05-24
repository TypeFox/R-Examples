iter1.norm <-
function(x,para0,lambda)
{
	sn=var(x)
	n=length(x)
	m=length(para0)/3
	xx=c()
	pdf.component=c()

	alpha=para0[1:m]
	mu=para0[(m+1):(2*m)]
	sigma=para0[(2*m+1):(3*m)]
	###E-step of EM-algorithm
	for (j in 1:m)
	{	
		pdf.component=cbind(pdf.component,dnorm(x,mu[j],sqrt(sigma[j]))*alpha[j]+1e-100/m)
		xx=cbind(xx,(x-mu[j])^2)
	}
	pdf=apply(pdf.component,1,sum)
	w=pdf.component/pdf
	###M-step of EM-algorithm
	alpha=(apply(w,2,sum)+lambda)/(n+m*lambda)
	mu=apply(w*x,2,sum)/apply(w,2,sum)
	sigma=(apply(w*xx,2,sum)+2*sn/sqrt(n))/(apply(w,2,sum)+2/n)
	
	###compute the log-likelihood and penalized log-likelihood value
	pdf.component=c()
	for (j in 1:m)
		pdf.component=cbind(pdf.component,dnorm(x,mu[j],sqrt(sigma[j]))*alpha[j]+1e-100/m)
	pdf=apply(pdf.component,1,sum)
	ln1=sum(log(pdf))
	pln1=ln1+sum(pn(sigma,sn,1/n))+lambda*sum(log(alpha))
	###output
	index=sort(mu,index.return = TRUE)$ix
	alpha=alpha[index]
	mu=mu[index]
	sigma=sigma[index]
	outpara=c(alpha,mu,sigma,ln1,pln1)
}
