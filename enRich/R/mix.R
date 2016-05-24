mix<-function(data, method=NULL, initialpara=NULL,fixoffset=FALSE, fixk=3,krange=c(0:10), exp.label=NULL, stopdiff=0.0001, parallel=TRUE)
{
## INPUT
## data: list of data, the first list is the region, which is a n x 3 matrix, where the name of the columns are  ("Chromosome", "Start", "Stop")
##          the second list contains the counts of ChIP experiments, which is a n x p matrix, where n is the number of regions and p is the number of experiments. 
##          At least one experiment count should be given.
## method: Can be "Poisson" or "NB" and it refers to the densities of the mixture distribution.
## initialpara: the initial parameters given for EM algorithm. In form of c("p", "lambda_S", "lambda_B") if method="Poisson" or c("p", "mu_S", "phi_S", "mu_B", "phi_B") if method="NB". 
##              Could be a matrix if initial values are the different for multiple experiments or a vector if initial values are the same. 
##              If not given, then a default value of (0.1, 10, 1) or (0.1, 10, 1, 1, 1) for method="Poisson" or "NB" respectively. 
## fixoffset: If TRUE, the offset of the signal distribution is fixed by the user and is the same for all experiments. 
##		  If FALSE, the offset is estimated empirically for each experiment.
## fixk: The value of the offset, when fixoffset = TRUE.

	data$count=as.matrix(data$count)
	data1=data$count
	Nexp=ncol(data1)
	if (is.null(exp.label)==1)
	{	
		if(is.null(colnames(data1))==1)
		{
			exp.label=paste("Experiment", 1:Nexp, sep="")
		}else
		{
			exp.label=colnames(data1)
		}
	}	
	colnames(data$count)=exp.label
	if (sum(method=="Poisson")+sum(method=="NB")==0)
	{
		stop('A method must be given, either "Poisson" or "NB"', call. = FALSE)
	}
	if (is.null(initialpara)==1)
	{
		if (method=="Poisson")
			initialpara=c(0.1, 10, 1)
		if (method=="NB")
			initialpara=c(0.1, 10, 1,1,1)
	}	
	if (method=="Poisson")
	{
		results<-matrix(0, Nexp, 4)
	}
	if (method=="NB")
	{
		results<-matrix(0, Nexp, 6)
	}
	if (Nexp==1)
	{
		cat("Only one experiment is being analysed.", '\n')
		if (fixoffset!=FALSE)
		{
			cat ("Fixed offset for each experiment equal to", fixk, '\n')
		}else
		{
			cat ("Offset is estimated by maximum likelihood in the range (", krange[1], ",", krange[length(krange)],")", '\n') 
		}
		results<-mix_single(data1, method=method, initialpara=initialpara, fixoffset=fixoffset, fixk=fixk, krange=krange,stopdiff=stopdiff)
		results=results[1:(length(results)-1)]
		results=t(as.matrix(results))
	}
	if (Nexp>1)
	{
		cat(Nexp, "experiments are modelled separately using latent mixture models",'\n')
		if (fixoffset!=FALSE)
		{
			cat ("Fixed offset for each experiment equal to", fixk, '\n')
		}else
		{
			cat ("Offset is estimated by maximum likelihood in the range (", krange[1], ",", krange[length(krange)],")", '\n') 
		}
		if (length(fixk)==1)
			fixk=rep(fixk, Nexp)
		if (is.vector(initialpara))
			initialpara=matrix(rep(initialpara, Nexp), Nexp, length(initialpara), byrow=T)
		jobs=paste("mix_single(data1[,",1:Nexp, "], method=method, initialpara=initialpara[", 1:Nexp, ",], fixoffset=fixoffset, fixk=fixk[", 1:Nexp, "], krange=krange, stopdiff=stopdiff)", sep="")
		if (parallel) 
		{		
			cl <- parallel::makeCluster(Nexp, type="PSOCK")			
			temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})
		}else
		{
			temp<-lapply(jobs, function(x) {eval(parse(text = x))})
		}
		for (i in 1:Nexp)
		{
			temp1=temp[[i]]
			results[i,]=temp1[1:(length(temp1)-1)]
		}	
	}
	if (parallel)
		parallel::stopCluster(cl)
	rownames(results)=exp.label
	if (method=="Poisson")
	{
		colnames(results)=c("p", "lambda_S", "lambda_B", "k")
	}
	if (method=="NB")
	{
		colnames(results)=c("p", "mu_S", "phi_S", "mu_B", "phi_B", "k")
	}
	object=list(data=data, parameters=results, method=method)	
	return(object)
}
