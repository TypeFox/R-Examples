separable1k<-
function (ymat, gamma = 1) 
{
	stopifnot(0==sum(is.na(as.vector(ymat))))
	n <- dim(ymat)[1]
	m <- dim(ymat)[2]
	o<-t(apply(ymat,1,sort))
	allmu<-matrix(NA,n,m-1)
	allsigma2<-matrix(NA,n,m-1)
	maxmu<-rep(-Inf,n)
	maxsig2<-rep(-Inf,n)
	
	for (j in 1:(m-1)){
		pr<-c(rep(1,j),rep(gamma,m-j))/(j+((m-j)*gamma))
		mu<-as.vector(o%*%pr)
		sigma2<-as.vector((o*o)%*%pr)-(mu*mu)
		chgmu<-(mu>maxmu)
		samemu<-(mu==maxmu)
		if (sum(chgmu)>0){
			maxmu[chgmu]<-mu[chgmu]
			maxsig2[chgmu]<-sigma2[chgmu]
		}
		if (sum(samemu)>0){
			maxsig2[samemu]<-pmax(sigma2[samemu],maxsig2[samemu])
		}
	}
	tstat <- as.vector(sum(ymat[,1]))
	expect <- sum(maxmu)
	vartotal <- sum(maxsig2)
	dev <- (tstat - expect)/sqrt(vartotal)
	pval <- 1 - pnorm(dev)
	list(pval = pval, deviate = dev, statistic = tstat, expectation = expect, 
		 variance = vartotal)
}
