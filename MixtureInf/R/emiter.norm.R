emiter.norm <-
function(x,para,beta,pens,para0,k)
{
	m=length(beta)
	n=length(x)
	sigma0=para0[[3]]
	alpha=para[1:(2*m)]
	mu=para[(2*m+1):(4*m)]
	sigma=para[(4*m+1):(6*m)]
	for (i in 1:k)
	{
		###E-step of EM-algorithm
		pdf.component=c()
		xx=c()
		for (j in 1:(2*m))
		{	
			pdf.component=cbind(pdf.component,dnorm(x,mu[j],sqrt(sigma[j]))*alpha[j]+1e-100/m)
			xx=cbind(xx,(x-mu[j])^2)
		}
		pdf=apply(pdf.component,1,sum)
		w=pdf.component/pdf

		###M-step of EM-algorithm
		ww=apply(w,2,sum)
		mu=apply(w*x,2,sum)/apply(w,2,sum)
		sigma=(apply(w*xx,2,sum)+2*pens[2]*sigma0)/(apply(w,2,sum)+2*pens[2])
		alpha=c()
		for (j in 1:m)
		{
			if (ww[2*j-1]/(ww[2*j-1]+ww[2*j])<=0.5)
				beta[j]=min((ww[2*j-1]+pens[1])/(ww[2*j-1]+ww[2*j]+pens[1]),0.5)
			else
				beta[j]=max((ww[2*j-1])/(ww[2*j-1]+ww[2*j]+pens[1]),0.5)
			alpha=c(alpha,(ww[2*j-1]+ww[2*j])/n*beta[j],(ww[2*j-1]+ww[2*j])/n*(1-beta[j]))
		}
	}
	###Compute the penalized log-likelihood value
	pdf.component=c()
	for (j in 1:(2*m))	
		pdf.component=cbind(pdf.component,dnorm(x,mu[j],sqrt(sigma[j]))*alpha[j]+1e-100/m)
	pdf=apply(pdf.component,1,sum)
	ln=sum(log(pdf))
	sigma0=array(rbind(sigma0,sigma0))
	pln=ln+pens[1]*sum(log(1-abs(1-2*beta)))+sum(pn(sigma,sigma0,pens[2]))
	output=c(alpha,mu,sigma,pln)
}
