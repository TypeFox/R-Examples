npmp <-
function(y, k=length(y), nrep, nb, alpha=1, theta=alpha, sigma=0, mixing_hyperprior= FALSE, basemeasure_hyperprior = FALSE,
	mixing_type="DP", algo="slice",  prior="gamma", a, b, a_a=1, b_a=1, 
	lb=NULL, ub=NULL, print = 1, ndisplay = nrep/4, plot.it = FALSE, pdfwrite = FALSE, ... )
{
	if(is.null(lb)) lb <- max(0,min(y)-10)
	if(is.null(ub)) ub <- max(y)+10
	if(prior=="gamma") priortype=1
	if(prior=="normal") priortype=2
	n <- length(y)
	if((mixing_type=="2PD") & (mixing_hyperprior==TRUE))
	{
		cat("Error: Random hyperparameters for the 2PD prior are not implemented\n")
		return()
	}
	if((algo=="polya-urn") & k<n)
	{
		cat("Warning: If algo='polya-urn', the upper bound for the number of components is fixed to the sample size n\n")
		k <- n
	}
	if(!((algo=="polya-urn") | (algo=="slice")))
	{
		cat("Error: Only 'slice' or 'polya-urn' are allowed values of 'algo'\n")
		return()
	}

	print <-  as.integer(table(factor(print, levels=1:5)))
	start.loop <- Sys.time()	
	res <- .C("npmpois",
				as.integer(c(n,k,ub+1,nrep, mixing_hyperprior, basemeasure_hyperprior, algo=="slice")),
				alpha=as.double(rep(alpha, nrep)),
				PDp_par=as.double(c(sigma, theta)),
				as.double(c(a,b,priortype,a_a,b_a)),
				as.double(y),
				as.integer(c(print,ndisplay)), 
				lambda=as.double(rep(1,k*nrep)), 
				pi=as.double(rep(c(1,rep(0,k-1)),nrep)),
				a =as.double(rep(a, nrep)),
				b =as.double(rep(b, nrep)),
				minslice=as.double(rep(1, nrep)),
				probability=as.double(rep(0,(ub+1)*nrep)),
				enne=as.double(rep(c(n,rep(0,k-1)),nrep)),
				PACKAGE="rmp"
				)
	cat("MCMC done\n")
	end.loop <- Sys.time()
	lambda<-matrix(res$lambda,byrow=T,nrow=nrep,ncol=k)
	pi<-matrix(res$pi,byrow=T,nrow=nrep,ncol=k)
	enne<-matrix(res$enne,byrow=T,nrow=nrep,ncol=k)	
	pmf<-matrix(res$probability,byrow=T,nrow=nrep,ncol=ub+1)
	
	post.pmf<-apply(pmf[(nb+1):nrep,],2,mean)
	sorted.p<-apply(pmf[(nb+1):nrep,],2,sort)
	pmf2.5 <-sorted.p[(0.025*(nrep-nb)),]
	pmf97.5 <-sorted.p[(0.975*(nrep-nb)),]
	
	sorted.n<-t(apply(enne[(nb+1):nrep,],1,sort,decreasing=T))
	n.mean <- apply(sorted.n,2,mean)
	arezero = function(x) sum(x==0)
	groups = k - mean( apply(sorted.n,1,arezero))

	post.lambda<-apply(lambda[(nb+1):nrep,],2,mean)
	post.pi<-apply(pi[(nb+1):nrep,],2,mean)
	post.alpha<-mean(res$alpha[(nb+1):nrep])

	empirical.pmf <- as.double(table(factor(y, levels=lb:ub))/length(y))
  
	if(plot.it){
	plot(lb:ub,empirical.pmf,type='h', main="", ylab="pmf", xlab="estimated pmf (red) and empirical pmf (black)") 
	points(lb:ub+.4,post.pmf[lb:ub+1],col=2,type='h')
	}
	
	if(pdfwrite){
		cat("\nWriting pdf files with traceplots ...\n")
		filename<-paste("traceplots ",date(),".pdf",sep="")
		pdf(filename)
		par(mfrow=c(3,3))
		plot(0,main ="TRACEPLOTS FOR LAMBDA -->")  
		for(i in 1:k) 	plot(lambda[,i],type='l') 
		plot(0,main ="TRACEPLOTS FOR PI -->")   
		for(i in 1:k)  plot(pi[,i],type='l') 
		par(mfrow=c(1,1))
		dev.off()

		filename<-paste("Estim pmf ",date(),".pdf",sep="")
		pdf(filename)
		plot(lb:ub,empirical.pmf,type='h',main="", ylab="pmf", xlab="estimated pmf (red) and empirical pmf (black)") 
		points(lb:ub+.4,post.pmf[lb:ub+1],col=2, type='h')
		dev.off()
		cat(paste("Pdf files with traceplots and posterior summaries written in ",getwd(),"\n\n",sep=""))
	}

    out<-structure(
		list(name = "DP mixture of Poisson kernels", mixing_type=mixing_type,
		  mcmc = list(nrep = nrep, nb = nb, time = (end.loop - start.loop)), 
		  mcmc.chains = list(alpha=res$alpha, lambda=lambda, pi=pi, a = res$a, b = res$b,
			minslice=res$minslice, size=enne, pmf=pmf),  
		  pmf = list(post.pmf=post.pmf[lb:ub+1], empirical=empirical.pmf, lower.95 = pmf2.5[lb:ub+1], 
		  upper.95 = pmf97.5[lb:ub+1], domain = lb:ub, lb=lb, ub=ub),
		  parameters = list(post.lambda=post.lambda, post.pi=post.pi, post.alpha=post.alpha),
  	            clustering = list(nmean=n.mean,groups=groups)), class  = "rmpobject")
	cat("\n")
    invisible(out)
    	}

dpmpoiss <- function(y, k=length(y), nrep, nb, alpha = 1, a, b, lb = NULL, ub = NULL, print = 1,
ndisplay = nrep/4, plot.it = FALSE, pdfwrite = FALSE, ...)
{
		npmp(y, k=k, nrep=nrep, nb=nb, alpha=alpha, sigma=0, mixing_hyperprior= FALSE, basemeasure_hyperprior = FALSE,
			mixing_type="DP", algo="slice",  prior="gamma", a=a, b=b, 
			lb=lb, ub=ub, print = print, ndisplay = ndisplay, plot.it = plot.it, pdfwrite = pdfwrite)
}