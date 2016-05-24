rmg <-
function(ydis, k=length(ydis), nrep, nb, alpha=1, theta=alpha, sigma=0, 
	mixing_hyperprior= FALSE, basemeasure_hyperprior = FALSE, mixing_type="DP", algo="slice",
	mu0=mean(ydis), kap=var(ydis), 
	atau=1, btau=2, a_a=1, b_a=1, lb=NULL, ub=NULL, print = 1, ndisplay = nrep/4,
	plot.it = FALSE, pdfwrite = FALSE, ...)
{
	
	n <- length(ydis)
	if((mixing_type=="2PD") & (mixing_hyperprior==TRUE))
	{
		cat("Error: Random hyperparameters for the 2PD prior are not implemented\n")
		return()
	}
	if((algo=="polya-urn") & k<n)
	{
		cat("Warning: If algo='polya-urn', the upper bound is fixed to the sample size n\n")
		k <- n
	}
	if(!((algo=="polya-urn") | (algo=="slice")))
	{
		cat("Error: Only 'slice' or 'polya-urn' are allowed values of 'algo'\n")
		return()
	}
	if(is.null(lb)) lb <- max(0,min(ydis)-10)
	if(is.null(ub)) ub <- max(ydis)+10
	print <-  as.integer(table(factor(print, levels=1:6)))
	start.loop <- Sys.time()			
	res <- .C("rmg"	,
			as.integer(c(n,k,ub+1,nrep,mixing_hyperprior,basemeasure_hyperprior,algo=="slice")),
			as.double(c(mu0, kap, atau, btau,a_a,b_a)),
			as.double(ydis),
			as.integer(c(print,ndisplay)), 
			alpha=as.double(rep(alpha, nrep)),
			PDp_par=as.double(c(sigma, theta)),
			mu0 = as.double(rep(mu0, nrep)),
			kappa = as.double(rep(kap, nrep)),
			atau = as.double(rep(atau, nrep)),
			btau = as.double(rep(btau, nrep)),
			mu=as.double(rep(mean(ydis),k*nrep)), 
			tau=as.double(rep(1/var(ydis),k*nrep)), 
			pi=as.double(rep(1/k,k*nrep)), 
			minslice=as.double(rep(1,nrep)), 
			probability=as.double(rep(0,(ub+1)*nrep)),
			enne=as.double(rep(1,k*nrep)),
			PACKAGE="rmp"
			)
	cat(nrep,"/",nrep, sep="")
	end.loop <- Sys.time()
	if(k==1)
	{
		mu<-res$mu
		tau<-res$tau
		pi<-res$pi
		post.mu<-mean(mu[(nb+1):nrep])
		post.tau<-mean(tau[(nb+1):nrep])
		post.pi<-1
		n.mean = n
		groups = 1
	}
	else
	{
		mu<-matrix(res$mu,byrow=T,nrow=nrep,ncol=k)
		tau<-matrix(res$tau,byrow=T,nrow=nrep,ncol=k)
		pi<-matrix(res$pi,byrow=T,nrow=nrep,ncol=k)
		enne<-matrix(res$enne,byrow=T,nrow=nrep,ncol=k)
		sorted.n<-t(apply(enne[(nb+1):nrep,],1,sort,decreasing=T))
		n.mean <- apply(sorted.n,2,mean)
		arezero = function(x) sum(x==0)
		groups = k - mean( apply(sorted.n,1,arezero))
		post.mu<-apply(mu[(nb+1):nrep,],2,mean)
		post.tau<-apply(tau[(nb+1):nrep,],2,mean)
		post.pi<-apply(pi[(nb+1):nrep,],2,mean)
	}

	cmf<-matrix(res$probability,byrow=T,nrow=nrep,ncol=ub+1)
	pmf<-cmf-cbind(rep(0,nrep),cmf[,1:(ncol(cmf)-1)])
	post.pmf<-try(apply(pmf[(nb+1):nrep,],2,mean))
	sorted.p<-try(apply(pmf[(nb+1):nrep,],2,sort))
	if(!is.null(dim(sorted.p)))
	{
	pmf2.5  <- apply(pmf[(nb+1):nrep,],2,quantile, prob=0.025)[lb:ub+1]
	pmf97.5 <- apply(pmf[(nb+1):nrep,],2,quantile, prob=0.975)[lb:ub+1]
	}
	else
	{
	pmf2.5 <- pmf97.5 <- NULL
	}
	empirical.pmf <- as.double(table(factor(ydis, levels=lb:ub))/length(ydis))

	post.alpha<-mean(res$alpha[(nb+1):nrep])

	if(plot.it)
	{
	plot(lb:ub,empirical.pmf,type='h', main="", ylab="pmf", xlab="estimated pmf (blue) and empirical pmf (black)") 
	points(lb:ub+.2,post.pmf[lb:ub+1], col=4, type='h')
	}
	
	if(pdfwrite)
	{
		cat("\nWriting pdf files with traceplots ...\n")
		filename<-paste("traceplots ",date(),".pdf",sep="")
		pdf(filename)
		plot(res$mu0,  type='l', xlab="MU_0")
		plot(res$kappa,type='l', xlab="KAPPA")		
		plot(res$atau, type='l', xlab="A_tau")		
		plot(res$btau, type='l', xlab="B_tau")		
		par(mfrow=c(3,3))
		if(k==1)
		{
		plot(mu,type='l', xlab="MU")
		plot(tau,type='l', xlab="TAU")
		}
		else
		{
		plot(0,main ="TRACEPLOTS FOR MU -->")  
		for(i in 1:k) 	plot(mu[,i],type='l')
		plot(0,main ="TRACEPLOTS FOR TAU -->")  
		for(i in 1:k)  plot(tau[,i],type='l')
		plot(0,main ="TRACEPLOTS FOR PI -->") 
		for(i in 1:k)  plot(pi[,i],type='l', ylim=c(0,1)) 
		plot(0,main ="TRACEPLOTS FOR SIZE -->") 
		for(i in 1:k)  plot(enne[,i],type='s', ylim=c(0,n)) 
		}
		par(mfrow=c(1,1))
		dev.off()

		if(plot.it)
		{
		filename<-paste("Estim pmf ",date(),".pdf",sep="")
		pdf(filename)
		plot(lb:ub,empirical.pmf,type='h', main="", ylab="pmf", xlab="estimated pmf (blue) and empirical pmf (black)") 
		points(lb:ub+.2,post.pmf[lb:ub+1],col=4,type='h')
		dev.off()
		}
		cat(paste("Pdf files with traceplots and posterior summaries written in ",getwd(),"\n\n",sep=""))
	}
    out<-structure(
		list(name = "Mixture of rounded Gaussians", mixing_type=mixing_type,
		mcmc = list(nrep = nrep, nb = nb, time = (end.loop - start.loop)), 
		mcmc.chains = list(alpha=res$alpha, mu0=res$mu0, kappa=res$kappa, 
		atau =res$atau, btau =res$btau, mu=mu, tau=tau, pi=pi, minslice=res$minslice, size=enne, cmf=cmf, pmf=pmf),
		pmf = list(post.pmf=post.pmf[lb:ub+1], empirical=empirical.pmf, 
		lower.95 = pmf2.5, upper.95 = pmf97.5, domain = lb:ub, lb=lb, ub=ub),
		parameters = list(post.mu=post.mu, post.tau=post.tau, post.pi=post.pi, post.alpha=post.alpha), 
  	          clustering = list(nmean=n.mean,groups=groups)), 
		class  = "rmpobject")
	cat("\n")
    invisible(out)
}
#
dpmrg <- function(ydis, k, nrep, nb, alpha, alpha_r = FALSE, mu0 = mean(ydis),
kap = var(ydis), atau, btau, a_a = 1, b_a = 1, lb = NULL, ub = NULL, print = 1,
ndisplay = nrep/4, plot.it = FALSE, pdfwrite = FALSE, ...)
{
	rmg(ydis, k=k, nrep, nb, alpha=alpha, theta=alpha, sigma=0, 
		mixing_hyperprior= alpha_r, basemeasure_hyperprior = FALSE,
		mu0=mu0, kap=kap, 
		atau=atau, btau=btau, a_a=a_a, b_a=b_a, lb=lb, ub=ub, print = print, ndisplay = ndisplay,
		plot.it = plot.it, pdfwrite = pdfwrite)
}
