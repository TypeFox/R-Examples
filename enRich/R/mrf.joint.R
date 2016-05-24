mrf.joint <-function(data, method=NULL, rep.vec=NULL, p.vec=NULL, exp.label=NULL, Niterations=10000, Nburnin=5000, Poisprior=NULL, NBprior=NULL, PoisNBprior=NULL, var.NB=NULL, var.q=NULL, parallel=TRUE)
{
## INPUT
## data: list of data, the first list is the region, which is a n x 3 matrix, where the name of the columns are ("chr", "start", "end")
##          the second list contains the counts of ChIP experiments, which is a n x p matrix, where n is the number of regions and p is the number of experiments. 
##          least one experiment count should be given.
## method: Can be "PoisNB", "Poisson" or "NB" and it refers to the densities of the mixture distribution.
## rep.vec: the vector of replicate indices, of length equal to the number of experiments. Technical replicates share the same index, 
##		e.g c(1,1, 2,2) for 4 experiments where the 1st and 2nd are two technical replicates and similarly the 3rd and 4th.
## p.vec: the vector of p indices, with p=P(X_s=1) for any region s. This vector is of length equal to the number of experiments. 
## exp.label: label of the experiments
## Niterations: the number of MCMC iteration steps. 
## Nburnin: the number of burn-in steps.
## Poisprior: the gamma priors for mean parameter lambda in Poisson-Poisson mixture, the first two are priors for signal and the second two are priors for background. 
##                       The prior should be a matrix 4 x p, each column represents the priors for each experiment.  Default values are (5, 1, 0.5, 1) for each single experiment. 
## NBprior: the gamma priors for mean mu and overdispersion parameters phi in NB-NB mixture, the first two are priors for mu_S for signal, the third and fourth are priors for phi_S 
##                    the fifth and sixth are priors for mu_B of background and the seventh and eighth are priors for phi_B. The prior should be
##                    a matrix 8 x p, each column represents the priors for each experiemnt. Default values are (5, 1, 1, 1, 0.5, 1, 1, 1) for each single experiment.
## PoisNBprior: the gamma priors for lambda_B and mu_S, phi_S in Poisson-NB mixture, the first two are priors for mu_S, the third and the fourth are priors for phi_S, 
##                      the fifth and the sixth are priors for lambda_B. The prior should be a matrix 6 x p, each column represents the priors for each  experiment.
##                      Default values are (5, 1, 1, 1, 0.5, 1) for each single experiment. 
## var.NB: the variances used in Metropolis-Hastling algorithm for estimates of (mu_S, phi_S, mu_B, phi_B) for NB mixture or for estimates of (mu_S, phi_S) for poisNB mixture. 
##              var.NB should be 4 x p or 2 x p matrix for NB and poisNB respectively, each column represents the variances used for each experiment.
##              Default values for each single experiment are (0.1, 0.1, 0.1, 0.1) or (0.1, 0.1) for NB and poisNB respectively. 
## var.q: the variances used in Metropolis-Hastling algorithm for estimates of q_0 and common ratio parameter when assume same p condition for multiple experiments. 
##           The number of components of var.q equals to number of experiment +1. Default values are 0.1 for each experiment and 0.3 for common ratio parameter. 
##           For example,  var.q=(0.1, 0.1, 0.3) for two experiments. 
	
	data$count=as.matrix(data$count)
	datacount=data$count
	Nexp=ncol(datacount)
	if (is.null(exp.label))
	{	
		if(is.null(colnames(datacount)))
		{
			exp.label=paste("Experiment", 1:Nexp, sep="")
		}else
		{
			exp.label=colnames(datacount)
		}
	}
	colnames(data$count)=exp.label
	if (sum(method=="Poisson")+sum(method=="NB")+sum(method=="PoisNB")==0)
	{
		stop('A method must be given, either "Poisson" or "NB" or "PoisNB" ' , call. = FALSE)
	}	
	if (length(rep.vec)>0 & length(rep.vec)!=Nexp)
	{
		stop("The length of rep.vec must be equal to the number of experiments", call. = FALSE)
	}	
	if (length(p.vec)>0 & length(p.vec)!=Nexp)
	{
		stop("The length of p.vec must be equal to the number of experiments", call. = FALSE)
	}
	if (any(p.vec==0)|any(rep.vec==0))
	{
		stop("Cannot use 0 in rep.vec or p.vec", call. = FALSE)
	}
	if (length(rep.vec)==0 & length(p.vec)==0)
	{
		stop("At least rep.vec or p.vec should be given", call. = FALSE)
	}
	if (length(rep.vec)>0&all(summary(factor(rep.vec))==1))## all indices of replicates are different
	{
		rep.vec=NULL
	}
	if (any(summary(factor(rep.vec))>1)& length(p.vec)>0)
	{	
		rlist=factor(rep.vec)
		rlevel=levels(rlist)
		rlevelvalue=summary(rlist)
		srlevel=subset(rlevel, rlevelvalue>1)
		srlevelvalue=subset(rlevelvalue, rlevelvalue>1)
		flist=factor(p.vec)
		level=levels(flist)
		levelvalue=summary(flist)
		splevel=subset(level, levelvalue>1)
		splevelvalue=subset(levelvalue, levelvalue>1)
		if ((Nexp-sum(srlevelvalue))<(Nexp-sum(splevelvalue)))
		{				
			stop ("Replicates should share the same p.vec index",call. = FALSE)
		}else if(length(rlevelvalue)==length(levelvalue)&&all(rlevelvalue==levelvalue))#all p.vec levels are the same as rep.vec levels
		{
			p.vec=NULL
		}
	}	
	if (is.null(Poisprior)==1)
	{
		Poisprior=c(5, 1, 0.5, 1)*matrix(1, 4, Nexp)
	}else if (is.vector(Poisprior))
	{
		Poisprior=as.matrix(Poisprior, 4, Nexp)
	}	
	if (is.null(NBprior)==1)
	{
		NBprior=c(5,1,1,1,0.5, 1,1,1)*matrix(1, 8, Nexp)
	}else if (is.vector(NBprior))
	{
		NBprior=as.matrix(NBprior, 8, Nexp)
	}
	if (is.null(PoisNBprior)==1)
	{
		PoisNBprior=c(5,1,1,1, 0.5,1)*matrix(1, 6, Nexp)
	}else if (is.vector(PoisNBprior))
	{
		PoisNBprior=as.matrix(PoisNBprior, 6, Nexp)
	}
	if (is.null(var.NB)==1)
	{
		if (method=="NB")
		{
			var.NB=c(0.1, 0.1, 0.1, 0.1)*matrix(1, 4, Nexp)
		}else
		{
			var.NB=c(0.1, 0.1)*matrix(1, 2, Nexp)
		}
	}else if (is.vector(var.NB))
	{
		var.NB=as.matrix(var.NB, length(var.NB)/Nexp, Nexp)
	}
	factor.chr=summary(factor(data$region[,1]))
	factor.name=names(factor.chr)
	Nchr=length(factor.chr)
	Ncl=Nchr
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
	paranumber=ifelse(method=="Poisson", 5, ifelse(method=="PoisNB", 6, 7))+1
	para.sample=list()
	para=vector("list", Nexp)
	names(para)=exp.label
	PP=matrix(0, length(datacount[,1]), Nexp)
	colnames(PP)=exp.label
	acrate=vector("list", Nexp)
	names(acrate)=exp.label
	tempdata=list()
	tempdata$region=data$region

##  separate modelling part p.vec and rep.vec are all different
	if ((length(p.vec)>0&all(summary(factor(p.vec))==1)))
	{	
		warning("All indices of p.vec are different: separate analyses for each experiment are conducted by calling mrf function", call. = FALSE)
		for (i in 1:Nexp)
		{
			cat(exp.label[i], "is treated as a single experiment, we use mrf function. \n")
			tempdata$count=data$count[,i]		
			tempresults=mrf(tempdata, method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior[,i], NBprior=NBprior[,i], PoisNBprior=PoisNBprior[,i], var.NB=var.NB[,i], parallel=parallel)
			para.sample[[i]]=tempresults$parameters.sample
			para[[i]]=tempresults$parameters
			PP[, i]=tempresults$PP	
			acrate[[i]]=tempresults$acrate.NB	
		}
		rep.vec=NULL
		p.vec=NULL
	}

## joint modelling part -- replicates
	if(!is.null(rep.vec)&is.null(p.vec))
	{
		rlist=factor(rep.vec)
		rlevel=levels(rlist)
		rlevelvalue=summary(rlist)
		srlevel=subset(rlevel, rlevelvalue>1)
		srlevelvalue=subset(rlevelvalue, rlevelvalue>1)
		for (i in 1:Nexp)
		{
			if (all(factor(rep.vec)[i]!=srlevel))
			{
				tempdata$count=data$count[,i]		
				cat (exp.label[i],"is treated as a single experiment, we use mrf function", '\n')
				tempresults=mrf(tempdata, method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior[,i], NBprior=NBprior[,i], PoisNBprior=PoisNBprior[,i], var.NB=var.NB[,i], parallel=parallel)
				para.sample[[i]]=tempresults$parameters.sample
				para[[i]]=tempresults$parameters
				PP[,i]=tempresults$PP	
				acrate[[i]]=tempresults$acrate.NB
			}	
		}
		for (j in 1: length(srlevelvalue))
		{
			Poisprior1=NULL
			PoisNBprior1=NULL
			NBprior1=NULL
			var.NB1=NULL
			data1=NULL
			exp.label1=NULL
			index=ifelse(factor(rep.vec)==srlevel[j], 1, 0)
			for (j1 in 1:Nexp)
			{
				if (index[j1]==1)
				{
					data1<-cbind(data1, datacount[,j1])
					NBprior1<-c(NBprior1, NBprior[,j1])
					Poisprior1<-c(Poisprior1, Poisprior[,j1])
					PoisNBprior1<-c(PoisNBprior1, PoisNBprior[,j1])
					var.NB1<-c(var.NB1, var.NB[,j1])
					exp.label1<-c(exp.label1, exp.label[j1])
				}
			}
			cat(exp.label1, "are treated as replicates and analyzed jointly. \n")
			datajoint=list()
			jobs=paste("mrf_rep(datajoint[[",1:Nchr, "]], method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior1, NBprior=NBprior1, PoisNBprior=PoisNBprior1, var.NB=var.NB1)", sep="")
			for (i in 1:Nchr)
			{
				datajoint[[i]]=subset(data1, data$region[,1]==names(factor.chr)[i])
			}
			if (parallel)
			{
				cl <- parallel::makeCluster(Ncl, type="PSOCK")			
				temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})			
			}else
			{
				temp<-lapply(jobs, function(x) {eval(parse(text = x))})
			}
			k1=0
			for (j1 in 1:Nexp)
			{					
				if (index[j1]==1)
				{
					k1=k1+1
					tempparameters=matrix(0, Nchr, paranumber)
					tempparameters.sample.temp=list()
					tempPP=NULL
					tempacrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
					for (i in 1:Nchr)
					{
						tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$parameters[[k1]], 2, mean), 4))
						tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$parameters[[k1]])
						tempPP=c(tempPP, temp[[i]]$PP)
						tempacrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
					}
					tempname=colnames(temp[[1]]$parameters[[1]])
					colnames(tempparameters)=c("Chromosome", tempname)
					tempparameters.sample=tempparameters.sample.temp[[1]]
					if(Nchr>2)
					{
						for (i in 2:Nchr)
						{
							tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
						}
					}										
					para.sample[[j1]]=tempparameters.sample
					para[[j1]]=tempparameters
					PP[,j1]=tempPP		
					acrate[[j1]]=tempacrate.NB			
				}
			}
		}
	}

## joint modelling part -- same p 
	if(is.null(rep.vec)&!is.null(p.vec)) 
	{
		flist=factor(p.vec)
		level=levels(flist)
		levelvalue=summary(flist)
		splevel=subset(level, levelvalue>1)
		splevelvalue=subset(levelvalue, levelvalue>1)
		for (i in 1:Nexp)
		{
			if (all(flist[i]!=splevel))
			{
				cat (exp.label[i],"is treated as a single experiment, we use mrf function", '\n')	
				tempdata$count=data$count[,i]				
				tempresults=mrf(tempdata, method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior[,i], NBprior=NBprior[,i], PoisNBprior=PoisNBprior[,i], var.NB=var.NB[,i], parallel=parallel)
				para.sample[[i]]=tempresults$parameters.sample
				para[[i]]=tempresults$parameters
				PP[,i]=tempresults$PP	
				acrate[[i]]=tempresults$acrate.NB			
			}
		}	
		for (j in 1: length(splevelvalue))
		{
			Poisprior1=NULL
			PoisNBprior1=NULL
			NBprior1=NULL
			var.NB1=NULL
			data1=NULL
			exp.label1=NULL
			index=ifelse(factor(p.vec)==splevel[j], 1, 0)
			for (j1 in 1:Nexp)
			{
				if (index[j1]==1)
				{
					data1<-cbind(data1, datacount[,j1])
					NBprior1<-c(NBprior1, NBprior[,j1])
					Poisprior1<-c(Poisprior1, Poisprior[,j1])
					PoisNBprior1<-c(PoisNBprior1, PoisNBprior[,j1])
					var.NB1<-c(var.NB1, var.NB[,j1])
					exp.label1<-c(exp.label1, exp.label[j1])
				}
			}
			cat(exp.label1, "are analyzed jointly where same p assumption are used. \n")
			datajoint=list()
			jobs=paste("mrf_sp(datajoint[[",1:Nchr, "]], method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior1, NBprior=NBprior1, PoisNBprior=PoisNBprior1, var.NB=var.NB1, var.q=var.q)", sep="")
			for (i in 1:Nchr)
			{
				datajoint[[i]]=subset(data1, data$region[,1]==names(factor.chr)[i])
			}
			if (parallel)
			{
				cl <- parallel::makeCluster(Ncl, type="PSOCK")			
				temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})			
			}else
			{
				temp<-lapply(jobs, function(x) {eval(parse(text = x))})
			}
			k1=0
			for (j1 in 1:Nexp)
			{					
				if (index[j1]==1)
				{
					k1=k1+1
					tempparameters=matrix(0, Nchr, paranumber)
					tempparameters.sample.temp=list()
					tempPP=NULL
					tempacrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
					for (i in 1:Nchr)
					{
						tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$parameters[[k1]], 2, mean), 4))
						tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$parameters[[k1]])
						tempPP=c(tempPP, temp[[i]]$PP[,k1])
						tempacrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
					}
					tempname=colnames(temp[[1]]$parameters[[1]])
					colnames(tempparameters)=c("Chromosome", tempname)
					tempparameters.sample=tempparameters.sample.temp[[1]]
					if (Nchr>1)
					{
						for (i in 2:Nchr)
						{
							tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
						}
					}										
					para.sample[[j1]]=tempparameters.sample
					para[[j1]]=tempparameters
					PP[,j1]=tempPP			
					acrate[[j1]]=tempacrate.NB				
				}
			}
		}
	}

## joint modelling part -- replicates and same p
	if (any(summary(factor(rep.vec))>1)& any(summary(factor(p.vec))>1))
	{	
		rlist=factor(rep.vec)
		rlevel=levels(rlist)
		rlevelvalue=summary(rlist)
		srlevel=subset(rlevel, rlevelvalue>1)
		srlevelvalue=subset(rlevelvalue, rlevelvalue>1)
		flist=factor(p.vec)
		level=levels(flist)
		levelvalue=summary(flist)
		splevel=subset(level, levelvalue>1)
		splevelvalue=subset(levelvalue, levelvalue>1)
		spindex=matrix(0, Nexp, length(splevelvalue))
		srindex=matrix(0, Nexp, length(srlevelvalue))
		if ((Nexp-sum(srlevelvalue))<(Nexp-sum(splevelvalue)))
			stop ("Replicates should share the same p.vec index",call. = FALSE)
		for (i in 1:Nexp)
		{
			if (all(flist[i]!=splevel))
			{
				tempdata$count=data$count[,i]
				cat (exp.label[i],"is treated as a single experiment, we use mrf function", '\n')	
				tempresults=mrf(tempdata, method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior[,i], NBprior=NBprior[,i], PoisNBprior=PoisNBprior[,i], var.NB=var.NB[,i], parallel=parallel)
				para.sample[[i]]=tempresults$parameters.sample
				para[[i]]=tempresults$parameters
				PP[,i]=tempresults$PP	
				acrate[[i]]=tempresults$acrate.NB					
			}
		}	
		if (length(srlevel)>0)
		{
			for (j in 1:length(srlevel))
			{
				srindex[,j]=ifelse(rlist==srlevel[j], 1, 0)
			}
		}
		for (j in 1:length(splevel))
		{
			spindex[,j]=ifelse(flist==splevel[j], 1, 0)
		}
		for (j in 1:length(splevel))
		{
			srdata=NULL
			srcode=NULL
			Nsr=NULL
			srNBprior=NULL
			srPoisprior=NULL
			srPoisNBprior=NULL
			srvar.NB=NULL
			srexp.label=NULL
			spexp.label=NULL
			finalindex=rep(0, Nexp)
			if (length(srlevel)>0)
			{
				for (k in 1:length(srlevel))
				{
					if(sum(srindex[,k]*spindex[,j])>1)
					{
						srdata=cbind(srdata, t(subset(t(datacount), srindex[,k]*spindex[,j]==1))) 
						srcode=c(srcode, subset(c(1:Nexp), srindex[,k]*spindex[,j]==1))
						srNBprior=cbind(srNBprior, t(subset(t(NBprior), srindex[,k]*spindex[,j]==1)))
						srPoisprior=cbind(srPoisprior, t(subset(t(Poisprior), srindex[,k]*spindex[,j]==1)))
						srPoisNBprior=cbind(srPoisNBprior, t(subset(t(PoisNBprior), srindex[,k]*spindex[,j]==1)))
						srvar.NB=cbind(srvar.NB, t(subset(t(var.NB), srindex[,k]*spindex[,j]==1)))
						srexp.label=c(srexp.label, subset(exp.label, srindex[,k]*spindex[,j]==1))
						Nsr=c(Nsr, sum(srindex[,k]*spindex[,j]))
						finalindex=finalindex+ifelse(srindex[,k]*spindex[,j]==1, 2, 0)
					}else
					{
						if(sum(subset(spindex[,j], srindex[,k]>0))>0)
						{
							stop ("Technical replicates should share the same p.vec index",call. = FALSE)
						}
					}
				}
			 }
			if (!is.null(srexp.label))
			{
				for (k in 1:length(srlevel))
				{
					tempcode=which(srindex[,k]==1)
					cat(exp.label[tempcode], "are treated as replicates and analyzed jointly. \n")
				}
				if (length(Nsr)>1)
				{
					cat(srexp.label, "are also analyzed jointly where same p assumption are used. \n")
				}
			}
			spdata=t(subset(t(datacount), spindex[,j]==1&rowSums(srindex)==0))
			spcode=subset(c(1:Nexp), spindex[,j]==1&rowSums(srindex)==0)
			spNBprior=t(subset(t(NBprior), spindex[,j]==1&rowSums(srindex)==0))
			spPoisprior=t(subset(t(Poisprior), spindex[,j]==1&rowSums(srindex)==0))
			spPoisNBprior=t(subset(t(PoisNBprior), spindex[,j]==1&rowSums(srindex)==0))
			spvar.NB=t(subset(t(var.NB), spindex[,j]==1&rowSums(srindex)==0))
			spexp.label=subset(exp.label, spindex[,j]==1&rowSums(srindex)==0)
			Nsp=sum(spindex[,j]==1&rowSums(srindex)==0)
			finalindex=finalindex+ifelse(spindex[,j]==1, 1, 0)
			spdata=as.matrix(spdata)
			if (!is.null(srexp.label)&&length(spexp.label)>0)
			{
				cat(srexp.label, "and", spexp.label, "are analyzed jointly where same p assumption are used. \n")
			}else if (length(spexp.label)>0)
			{
				cat(spexp.label, "are analyzed jointly where same p assumption are used. \n")
			}
## 		
			if (length(Nsr)==0)
			{
				nsr1=0
				nsr2=0
				nsp1=0
				nsp2=Nsp
				data1=NULL
				Poisprior1=NULL
				NBprior1=NULL
				PoisNBprior1=NULL
				exp.label1=NULL
				var.NB1=NULL
				data2 =spdata
				Poisprior2=as.vector(spPoisprior)
				NBprior2=as.vector(spNBprior)
				PoisNBprior2=as.vector(spPoisNBprior)
				exp.label2=spexp.label		
				var.NB2=as.vector(spvar.NB)						
			}  	
			if (length(Nsr)==1)
			{
				nsr1=Nsr
				nsr2=0
				nsp2=Nsp
				nsp1=0
				data1=srdata
				Poisprior1=as.vector(srPoisprior)
				NBprior1=as.vector(srNBprior)
				PoisNBprior1=as.vector(srPoisNBprior)
				exp.label1=srexp.label
				var.NB1=as.vector(srvar.NB)
				if (nsp2>0)
				{
					data2 =spdata
					Poisprior2=as.vector(spPoisprior)
					NBprior2=as.vector(spNBprior)
					PoisNBprior2=as.vector(spPoisNBprior)
					exp.label2=spexp.label
					var.NB2=as.vector(spvar.NB)
				}else
				{
					data2=NULL
					Poisprior2=NULL
					PoisNBprior2=NULL
					NBprior2=NULL
					exp.label2=NULL
				}				
			}  			
			if (length(Nsr)==2)
			{
				nsr1=Nsr[1]
				nsr2=Nsr[2]
				nsp2=Nsp
				nsp1=0
				data1=srdata[, 1:nsr1]
				Poisprior1=as.vector(srPoisprior[,1:nsr1])
				NBprior1=as.vector(srNBprior[,1:nsr1])
				PoisNBprior1=as.vector(srPoisNBprior[,1:nsr1])
				var.NB1=as.vector(srvar.NB[,1:nsr1])
				exp.label1=srexp.label[1:nsr1]
				if (nsp2>0)
				{
					data2=cbind(srdata[, (nsr1+1):(nsr1+nsr2)], spdata)
					Poisprior2=c(as.vector(srPoisprior[, (nsr1+1):(nsr1+nsr2)]), spPoisprior)
					NBprior2=c(as.vector(srNBprior[ ,(nsr1+1):(nsr1+nsr2)]), spPoisprior)
					PoisNBprior2=c(as.vector(srPoisNBprior[, (nsr1+1):(nsr1+nsr2)]), spPoisNBprior)
					var.NB2=c(as.vector(srvar.NB[, (nsr1+1):(nsr1+nsr2)]), spvar.NB)
					exp.label2=c(srexp.label[(nsr1+1):(nsr1+nsr2)], spexp.label)
				}else
				{
					data2=srdata[, (nsr1+1):(nsr1+nsr2)]
					Poisprior2=as.vector(srPoisprior[, (nsr1+1):(nsr1+nsr2)])
					NBprior2=as.vector(srNBprior[ ,(nsr1+1):(nsr1+nsr2)])
					PoisNBprior2=as.vector(srPoisNBprior[, (nsr1+1):(nsr1+nsr2)])
					var.NB2=as.vector(srvar.NB[, (nsr1+1):(nsr1+nsr2)])
					exp.label2=srexp.label[(nsr1+1):(nsr1+nsr2)]
				}
			}  
			exp.labelf=c(exp.label1, exp.label2)

##
			if (is.null(data2)) ## only use replicates
			{
				datajoint=list()
				jobs=paste("mrf_rep(datajoint[[",1:Nchr, "]], method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior1, NBprior=NBprior1, PoisNBprior=PoisNBprior1, var.NB=var.NB1)", sep="")
				for (i in 1:Nchr)
				{
					datajoint[[i]]=subset(data1, data$region[,1]==names(factor.chr)[i])
				}
				if (parallel)
				{
					cl <- parallel::makeCluster(Ncl, type="PSOCK")			
					temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})			
				}else
				{
					temp<-lapply(jobs, function(x) {eval(parse(text = x))})
				}
				for (j1 in 1:nsr1)
				{					
					tempparameters=matrix(0, Nchr, paranumber)
					tempparameters.sample.temp=list()
					tempPP=NULL
					tempacrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
					for (i in 1:Nchr)
					{
						tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$parameters[[j1]], 2, mean), 4))
						tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$parameters[[j1]])
						tempPP=c(tempPP, temp[[i]]$PP)
						tempacrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
					}
					tempname=colnames(temp[[1]]$parameters[[1]])
					colnames(tempparameters)=c("Chromosome", tempname)
					tempparameters.sample=tempparameters.sample.temp[[1]]
					if (Nchr>1)
					{
						for (i in 2:Nchr)
						{
							tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
						}
					}										
					para.sample[[srcode[j1]]]=tempparameters.sample
					para[[srcode[j1]]]=tempparameters
					PP[,srcode[j1]]=tempPP		
					acrate[srcode[j1]]=tempacrate.NB					
				}
			}
				
			if (is.null(data1)) ## only use same p
			{
				datajoint=list()
				jobs=paste("mrf_sp(datajoint[[",1:Nchr, "]], method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior=Poisprior1, NBprior=NBprior1, PoisNBprior=PoisNBprior1, var.NB=var.NB1, var.q=var.q)", sep="")
				for (i in 1:Nchr)
				{
					datajoint[[i]]=subset(data2, data$region[,1]==names(factor.chr)[i])
				}
				if (parallel)
				{
					cl <- parallel::makeCluster(Ncl, type="PSOCK")			
					temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})			
				}else
				{
					temp<-lapply(jobs, function(x) {eval(parse(text = x))})
				}
				for (j1 in 1:nsp2)
				{					
					tempparameters=matrix(0, Nchr, paranumber)
					tempparameters.sample.temp=list()
					tempPP=NULL
					tempacrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
					for (i in 1:Nchr)
					{
						tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$parameters[[j1]], 2, mean), 4))
						tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$parameters[[j1]])
						tempPP=c(tempPP, temp[[i]]$PP[,j1])
						tempacrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
					}
					tempname=colnames(temp[[1]]$parameters[[1]])
					colnames(tempparameters)=c("Chromosome", tempname)
					tempparameters.sample=tempparameters.sample.temp[[1]]
					if (Nchr>1)
					{
						for (i in 2:Nchr)
						{
							tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
						}
					}										
					para.sample[[spcode[j1]]]=tempparameters.sample
					para[[spcode[j1]]]=tempparameters
					PP[,spcode[j1]]=tempPP		
					acrate[spcode[j1]]=tempacrate.NB					
				}
			}

			if (!is.null(data1)&&!is.null(data2)) ## using both case
			{
				datajoint1=list()
				datajoint2=list()
				jobs=paste("mrf_srsp(datajoint1[[",1:Nchr, "]], datajoint2[[",1:Nchr, "]], nsr1=nsr1, nsr2=nsr2, nsp1=nsp1, nsp2=nsp2, method=method, Niterations=Niterations, Nburnin=Nburnin, Poisprior1=Poisprior1, NBprior1=NBprior1, PoisNBprior1=NBprior1, Poisprior2=Poisprior2, NBprior2=NBprior2, PoisNBprior2=PoisNBprior2, var.NB1=var.NB1, var.NB2=var.NB2)", sep="")
				for (i in 1:Nchr)
				{
					datajoint1[[i]]=subset(data1, data$region[,1]==names(factor.chr)[i])
					datajoint2[[i]]=subset(data2, data$region[,1]==names(factor.chr)[i])
				}
				if (parallel)
				{
					cl <- parallel::makeCluster(Ncl, type="PSOCK")			
					temp<-parallel::clusterApplyLB(cl, jobs, function(x) {eval(parse(text = x))})			
				}else
				{
					temp<-lapply(jobs, function(x) {eval(parse(text = x))})
				}
				tempacrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
				for (i in 1:Nchr)
				{
					tempacrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
				}
				if (nsr1>0)
				{
					for (j1 in 1:nsr1)
					{					
						tempparameters=matrix(0, Nchr, paranumber)
						tempparameters.sample.temp=list()
						tempPP=NULL
						for (i in 1:Nchr)
						{
							tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$para.con1[[j1]], 2, mean), 4))
							tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$para.con1[[j1]])
							tempPP=c(tempPP, temp[[i]]$PP.con1[,j1])
						}
						tempname=colnames(temp[[1]]$para.con1[[1]])
						colnames(tempparameters)=c("Chromosome", tempname)
						tempparameters.sample=tempparameters.sample.temp[[1]]
						if (Nchr>1)
						{
							for (i in 2:Nchr)
							{
								tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
							}
						}										
						para.sample[[srcode[j1]]]=tempparameters.sample
						para[[srcode[j1]]]=tempparameters
						PP[,srcode[j1]]=tempPP	
						acrate[srcode[j1]]=tempacrate.NB	
					}
				}
				if (nsp1>0)
				{
					for (j1 in 1:nsp1)
					{					
						tempparameters=matrix(0, Nchr, paranumber)
						tempparameters.sample.temp=list()
						tempPP=NULL						
						tempacrate.NB=matrix(0, Nchr, (length(temp[[1]]$acrate.NB)+1))
						for (i in 1:Nchr)
						{
							tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$para.con1[[j1+nsr1]], 2, mean), 4))
							tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$para.con1[[j1+nsr1]])
							tempPP=c(tempPP, temp[[i]]$PP.con1[,j1+nsr1])
							tempacrate.NB[i,]=c(factor.name[i], temp[[i]]$acrate.NB)
						}
						tempname=colnames(temp[[1]]$para.con1[[1]])
						colnames(tempparameters)=c("Chromosome", tempname)
						tempparameters.sample=tempparameters.sample.temp[[1]]
						if (Nchr>1)
						{
							for (i in 2:Nchr)
							{
								tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
							}
						}										
						para.sample[[spcode[j1]]]=tempparameters.sample
						para[[spcode[j1]]]=tempparameters
						PP[,spcode[j1]]=tempPP		
						acrate[spcode[j1]]=tempacrate.NB					
					}
				}
				if (nsr2>0)
				{
					for (j1 in 1:nsr2)
					{					
						tempparameters=matrix(0, Nchr, paranumber)
						tempparameters.sample.temp=list()
						tempPP=NULL						
						for (i in 1:Nchr)
						{
							tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$para.con2[[j1]], 2, mean), 4))
							tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$para.con2[[j1]])
							tempPP=c(tempPP, temp[[i]]$PP.con2[,j1])
						}
						tempname=colnames(temp[[1]]$para.con2[[1]])
						colnames(tempparameters)=c("Chromosome", tempname)
						tempparameters.sample=tempparameters.sample.temp[[1]]
						if (Nchr>1)
						{
							for (i in 2:Nchr)
							{
								tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
							}
						}										
						para.sample[[srcode[j1+nsr1]]]=tempparameters.sample
						para[[srcode[j1+nsr1]]]=tempparameters
						PP[,srcode[j1+nsr1]]=tempPP		
						acrate[srcode[j1+nsr1]]=tempacrate.NB					
					}
				}
				if (nsp2>0)
				{
					for (j1 in 1:nsp2)
					{					
						tempparameters=matrix(0, Nchr, paranumber)
						tempparameters.sample.temp=list()
						tempPP=NULL
						for (i in 1:Nchr)
						{
							tempparameters[i,]=c(factor.name[i], round(apply(temp[[i]]$para.con2[[j1+nsr2]], 2, mean), 4))
							tempparameters.sample.temp[[i]]=cbind(chr=factor.name[i], temp[[i]]$para.con2[[j1+nsr2]])
							tempPP=c(tempPP, temp[[i]]$PP.con2[,j1+nsr2])
						}
						tempname=colnames(temp[[1]]$para.con2[[1]])
						colnames(tempparameters)=c("Chromosome", tempname)
						tempparameters.sample=tempparameters.sample.temp[[1]]
						for (i in 2:Nchr)
						{
							tempparameters.sample=rbind(tempparameters.sample, tempparameters.sample.temp[[i]])
						}										
						para.sample[[spcode[j1+nsp1]]]=tempparameters.sample
						para[[spcode[j1+nsp1]]]=tempparameters
						PP[,spcode[j1+nsp1]]=tempPP		
						acrate[spcode[j1+nsp1]]=tempacrate.NB					
					}
				}
			}
		}
	}	

## Main OUTPUT
##  parameters: the list of mean of sample matrix of parameters for experiments. 
##  parameters.sample: the list of parameters for experiments and each list is the samples matrix drawing from the posterior distributions of parameters. 
##                      The samples are collected one from every ten steps right after burn-in step. 
##                     The column names for the matrix are (q_1, q_0, lambda_S, pi, lambda_B) if method="Poisson" or (q_1, q_0, mu_S, phi_S, pi, mu_B, phi_B) if method ="NB" 
##                      or (q_1, q_0, mu_S, phi_S, pi, lambda_B) if method="PoisNB", where q_1 and q_0 are the probability of  current region is enrich given the previous region
##                      is enrich or not enrich, respectively. lambda_S is mean of Poisson distribution for signal; mu_S and phi_S are mean and overdispersion of NB distribution for signal;
##                      lambda_B is mean of Poisson distribution for background; mu_B and phi_B are mean and overdispersion of NB distribution for background; 
##                      pi is the zero portion in zero-inflated model for background. 
## PP: the matrix of posterior probability for experiments and each column contains the posterior probability that the n region are enriched. 
	result=list(data=data, parameters=para, parameters.sample=para.sample, PP=PP, rep.vec=rep.vec, p.vec=p.vec, method=method, acrate.NB=acrate)
	return(result)
}
