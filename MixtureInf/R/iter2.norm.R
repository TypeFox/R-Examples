iter2.norm <-
function(x,para,beta,an,para0)
{
	m=length(beta)
	mu0=para0[[2]]
	sigma0=para0[[3]]
	alpha=para[1:(2*m)]
	mu=para[(2*m+1):(4*m)]
	sigma=para[(4*m+1):(6*m)]
	pdf.component=c()
	xx=c()

	###Calculate the cut points of parameter space of mu
	eta=rep(0,(m+1))
	eta[1]=min(x)
	eta[m+1]=max(x)
	if (m>1)
		eta[2:m]=(mu0[1:(m-1)]+mu0[2:m])/2
	
	###E-step of EM-algorithm
	for (j in 1:(2*m))
	{	
		pdf.component=cbind(pdf.component,dnorm(x,mu[j],sqrt(sigma[j]))*alpha[j]+1e-100/m)
		xx=cbind(xx,(x-mu[j])^2)
	}
	pdf=apply(pdf.component,1,sum)
	w=pdf.component/pdf
	
	###M-step of EM-algorithm
	alpha=apply(w,2,mean)
	mu=apply(w*x,2,sum)/apply(w,2,sum)
	sigma=(apply(w*xx,2,sum)+2*an*sigma0)/(apply(w,2,sum)+2*an)
	alpha1=c()
	mu1=c()
	for (j in 1:m)
	{
	alpha1=c(alpha1,(alpha[2*j-1]+alpha[2*j])*beta[j],(alpha[2*j-1]+alpha[2*j])*(1-beta[j]))
	mu1=c(mu1,min(max(mu[2*j-1],eta[j]),eta[j+1]),min(max(mu[2*j],eta[j]),eta[j+1]))
	}
	
	###Compute the log-likelihood value
	pdf.component=c()
	for (j in 1:(2*m))	
		pdf.component=cbind(pdf.component,dnorm(x,mu1[j],sqrt(sigma[j]))*alpha1[j]+1e-100/m)
	pdf=apply(pdf.component,1,sum)	
	ln=sum(log(pdf))
	sigma0=array(rbind(sigma0,sigma0))
	pln=ln+sum(pn(sigma,sigma0,an))
	output=c(alpha1,mu1,sigma,pln)
}
