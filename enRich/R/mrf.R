mrf<-function(data, method=NULL, exp.label=NULL, Niterations=10000, Nburnin=5000, Poisprior=c(5, 1, 0.5, 1), NBprior=c(5, 1, 1, 1, 0.5, 1, 1, 1), PoisNBprior=c(5,1,1,1, 0.5,1), var.NB=c(0.1, 0.1, 0.1, 0.1), parallel=TRUE)
{	
## INPUT
## data: list of data, the first list is the region, which is a n x 3 matrix, where the name of the columns are ("chr", "start", "end")
##          the second list contains the counts of a single ChIP experiment, which is a vector of size n. where n is the number of regions and p is the number of experiments. 
## method: Can be "PoisNB", "Poisson" or "NB" and it refers to the densities of the mixture distribution.
## exp.label: label of the experiments
## Niterations: the number of MCMC iteration steps. 
## Nburnin: the number of burn-in steps.
## Poisprior: the gamma priors for mean parameter lambda in Poisson-Poisson mixture, the first two are priors for signal and the second two are priors for background. 
##                 Default values are (5,1, 0.5, 1). 
## NBprior: the gamma priors for mean mu and overdispersion parameters phi in NB-NB mixture, the first two are priors for mu_S for signal, the third and fourth are priors for phi_S 
##                    the fifth and sixth are priors for mu_B of background and the seventh and eighth are priors for phi_B. Default values are (5, 1, 1, 1, 0.5, 1, 1, 1).
## PoisNBprior: the gamma priors for lambda_B and mu_S, phi_S in Poisson-NB mixture, the first two are priors for mu_S, the third and the fourth are priors for phi_S, 
##                      the fifth and the sixth are priors for lambda_B. Default values are (5, 1,1,1, 0.5, 1). 
## var.NB: the variances used in Metropolis-Hastling algorithm for estimates of (mu_S, phi_S, mu_B, phi_B) for NB mixture or for estimates of (mu_S, phi_S) for poisNB mixture. 
##              Default values are (0.1, 0.1, 0.1, 0.1) or (0.1, 0.1) for NB and poisNB respectively. 

	data$count=as.matrix(data$count)
	if (ncol(data$count)>1)
	{
		stop("mrf is for single ChIP experiment, use mrf.joint for more than one ChIP experiments")
	}
	data2=data$count
	if (sum(method=="Poisson")+sum(method=="NB")+sum(method=="PoisNB")==0)
	{
		stop('A method must be given, either "Poisson" or "NB" or "PoisNB" ' , call. = FALSE)
	}
	if (is.null(exp.label))
	{	
		if(is.null(colnames(data2)))
		{
			exp.label="Experiment1"
		}else
		{
			exp.label=colnames(data1)		
		}
	}
	colnames(data$count)=exp.label
	if (Niterations<=Nburnin)
		stop("Niteration must larger than Nburnin", call.=FALSE)
	if (method=="PoisNB")
	{
		met=0
		if (length(var.NB)>=2)
		{
			var.NB=c(var.NB[1], var.NB[2], 0, 0)
		}
		Poisprior=c(0,0, PoisNBprior[5:6])
		NBprior=c(PoisNBprior[1:4], 0,0,0,0)
	}
	if (method=="Poisson")
	{
		met=1
		NBprior=rep(0, 8)
	}
	if (method=="NB")
	{
		met=2
		Poisprior=rep(0, 4)
	}
	paranumber=ifelse(method=="Poisson", 5, ifelse(method=="PoisNB", 6, 7))+1
	factor.chr=summary(factor(data$region[,1]))
	factor.name=names(factor.chr)
	Nchr=length(factor.chr)
	if (Nchr>1&parallel)
	{
		parallel=TRUE
	}else
	{
		parallel=FALSE
	}
	if (parallel && !("parallel" %in% names(sessionInfo()$otherPkgs))) 
	{
		message("\nPackage 'parallel' is not loaded - parallel processing disabled. Please load multicore with library(parallel). \n")
		parallel = FALSE
	}
	data1=list()
	jobs=paste("mrf_single(data1[[",1:Nchr, "]], met=met, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior, NBprior=NBprior, PoisNBprior=PoisNBprior, var.NB=var.NB)", sep="")
	for (i in 1:Nchr)
	{
		data1[[i]]=subset(data2, data$region[,1]==names(factor.chr)[i])
	}
	if (parallel)
	{
		cl <- parallel::makeCluster(Nchr, type="PSOCK")			
		temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})
	}else
	{
		temp<-lapply(jobs, function(x) {eval(parse(text = x))})
	}
	if (parallel)
		parallel::stopCluster(cl)

## OUTPUT
## parameters: the mean of parameters sample for experiment. 
## parameters.sample: the samples matrix drawing from the posterior distributions of parameters. The samples are collected one from every ten steps right after burn-in step.         
##                    The column names for the matrix are (q_1, q_0, lambda_S, pi, lambda_B) if method="Poisson" or (q_1, q_0, mu_S, phi_S, pi, mu_B, phi_B) if method ="NB" 
##                    or (q_1, q_0, mu_S, phi_S, pi, lambda_B) if method="PoisNB", where q_1 and q_0 are the probability of  current region is enrich given the previous region
##                    is enrich or not enrich, respectively. lambda_S is mean of Poisson distribution for signal; mu_S and phi_S are mean and overdispersion of NB distribution for signal;
##                    lambda_B is mean of Poisson distribution for background; mu_B and phi_B are mean and overdispersion of NB distribution for background; 
##                    pi is the zero portion in zero-inflated model for background. 
## PP: the list of posterior probability for experiments and each list contains the posterior probability that the n region are enriched. 
	parameters=matrix(0, Nchr, paranumber)
	parameters.sample.temp=list()
	PP=NULL
	acrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
	for (i in 1:Nchr)
	{
		parameters[i,]=c(factor.name[i], temp[[i]]$parameters)
		parameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$parameters.sample)
		PP=c(PP, temp[[i]]$PP)
		acrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
	}
	temp.parameters.name=colnames(temp[[1]]$parameters)
	colnames(parameters)=c("Chromosome", temp.parameters.name)
	parameters.sample=parameters.sample.temp[[1]]
	if (Nchr>1)
	{
		for (i in 2:Nchr)
		{
			parameters.sample=rbind(parameters.sample, parameters.sample.temp[[i]])
		}
	}
	results=list(data=data, parameters=parameters, parameters.sample=parameters.sample, PP=PP, method=method, acrate.NB=acrate.NB)
	return(results)
}
