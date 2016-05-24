mix.joint <-function(data, method=NULL, para.sep=NULL, rep.vec=NULL, p.vec=NULL, exp.label=NULL, stopdiff=0.0001)
{
## INPUT
## data: list of data, the first list is the region, which is a n x 3 matrix, where the name of the columns are  ("Chromosome", "Start", "Stop")
##          the second list contains the counts of ChIP experiments, which is a n x p matrix, where n is the number of regions and p is the number of experiments. 
##          At least two experiments count should be given.
## para.sep: the initial parameter used for mix.joint fitting. It could be the result of the mix function or given by user. 
## method: Can be "Poisson" or "NB" and it refers to the densities of the mixture distribution.
## rep.vec: the vector of replicate indices, of length equal to the number of experiments. Technical replicates share the same index, 
##		e.g c(1,2,2,3,4,4,5,6) for 8 experiments where the 2nd and 3rd are two technical replicates and similarly the 5th and 6th.
## p.vec: the vector of p indices, with p=P(X_s=1) for any region s. This vector is of length equal to the number of experiments. 
##	    Experiments with the same probability of enrichement share the same p index, e.g. technical replicates and/or proteins with
##	    a similar number of binding sites, e.g. c(1,1,1,2,3,3,3,4) if the first three experiments have the same p and similarly the 5th, 6th and 7th experiments.
##        This allows to propertly account for the different IP efficiencies in the joint analysis. 
## 	    At least one of rep.vec or p.vec should be given. For those experiments which do not share the same index (p.vec or rep.vec) with any other experiments a single mixture model will be fitted. 

	data$count=as.matrix(data$count)
	Nexp=ncol(data$count)
	if (Nexp<2)
	{
		stop("There is only one experiment, cannot do a joint analysis!", call. = FALSE)
	}
	datacount=data$count
	if (is.null(exp.label)==1)
	{	
		if(is.null(colnames(data$count))==1)
		{
			exp.label=paste("Experiment", 1:Nexp, sep="")
		}else
		{
			exp.label=colnames(data$count)
		}
	}
	colnames(data$count)=exp.label
	if (sum(method=="Poisson")+sum(method=="NB")==0)
	{
		stop('A method must be given, either "Poisson" or "NB"', call. = FALSE)
	}	
	if (is.null(para.sep)==1)
	{
		stop('No values are given in para.sep, need to run mix first.', call. = FALSE)
	}
	initialpara=para.sep
	results=para.sep
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
			stop ("Replicates should share the same p.vec index",call. = FALSE)
	}
	if (length(p.vec)>0&all(summary(factor(p.vec))==1))
	{	
		warning("All indices of p.vec are different: separate analyses for each experiment are conducted", call. = FALSE)
		results=initialpara
	}
	if (length(p.vec)==0)
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
				cat ("Experiment", i,"is treated as a single experiment", '\n')				
			}
		}
		for (j in 1: length(srlevelvalue))
		{
			initialpara1=NULL
			data1=NULL
			index=ifelse(factor(rep.vec)==srlevel[j], 1, 0)
			for (j1 in 1:Nexp)
			{
				if (index[j1]==1)
				{
					initialpara1=c(initialpara1, initialpara[j1, ])
					data1<-cbind(data1, datacount[,j1])	
				}
			}
			tempresults=mix_srsp1(datasr=data1, datasp=NULL, method=method, initialpara_sr=initialpara1, initialpara_sp=NULL, Nsr=levelvalue[j], Nsp=0)$para_sr
			k=0
			for (j1 in 1:length(srlevelvalue))
			{
				if (index[j1]==1)
				{
					k=k+1
					results[j1,]=tempresults[k,]
					
				}
			}
		}
	}
	if (length(rep.vec)==0& any(summary(factor(p.vec))>1))
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
				cat ("Experiment", i,"is treated as a single experiment", '\n')				
			}
		}	
		for (j in 1:length(splevel))
		{
			initialpara1=NULL
			data1=NULL
			index=ifelse(flist==splevel[j], 1, 0)
			for (j1 in 1:Nexp)
			{
				if (index[j1]==1)
				{
					initialpara1=c(initialpara1, initialpara[j1, ])
					data1<-cbind(data1, datacount[,j1])	
				}
			}
			tempresults=mix_srsp1(datasr=NULL, datasp=data1, method=method, initialpara_sr=NULL, initialpara_sp=initialpara1, Nsr=0, Nsp=levelvalue[j])$para_sp
			k=1
			for (j1 in 1:Nexp)
			{
				if (index[j1]==1)
				{
					results[j1,]=tempresults[k,]
					k=k+1
				}
			}
		}
	}
	if (length(rep.vec)>0& any(summary(factor(p.vec))>1))
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
				cat ("Experiment", i,"is treated as a single experiment", '\n')	
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
			srinitialpara=NULL
			Nsr=NULL
			finalindex=rep(0, Nexp)
			if (length(srlevel)>0)
			{
				for (k in 1:length(srlevel))
				{
					if(sum(srindex[,k]*spindex[,j])>1)
					{
						srdata=cbind(srdata, t(subset(t(datacount), srindex[,k]*spindex[,j]==1))) 
						srinitialpara=rbind(srinitialpara, subset(initialpara, srindex[,k]*spindex[,j]==1))
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
			spdata=t(subset(t(datacount), spindex[,j]==1&rowSums(srindex)==0))
			spinitialpara=subset(initialpara,spindex[,j]==1&rowSums(srindex)==0)
			Nsp=sum(spindex[,j]==1&rowSums(srindex)==0)
			finalindex=finalindex+ifelse(spindex[,j]==1, 1, 0)
			srinitialpara1=NULL
			spinitialpara1=NULL
			if (length(srinitialpara)>0)
			{
				for (k in 1:nrow(srinitialpara))
				{
					srinitialpara1=c(srinitialpara1, srinitialpara[k,])
				}
			}else
			{
				Nsr=0
			}
			if (length(spinitialpara)>0)
			{
				for (k in 1:nrow(spinitialpara))
				{
					spinitialpara1=c(spinitialpara1, spinitialpara[k,])
				}
			}else
			{
				Nsp=0
			}
			tempresults=mix_srsp1(datasr=srdata, datasp=spdata, initialpara_sr=srinitialpara1, initialpara_sp=spinitialpara1, method=method, Nsr=Nsr, Nsp=Nsp)
			tempsrresults=tempresults$para_sr
			tempspresults=tempresults$para_sp
			srk=0
			spk=0
			for (j1 in 1:Nexp)
			{
				if (finalindex[j1]==3)
				{
					srk=srk+1
					results[j1,]=tempsrresults[srk,]
				}else
				{
					if (finalindex[j1]==1)
					{
						spk=spk+1
						results[j1,]=tempspresults[spk,]
					}
				}
			}
		}
	}
	rownames(results)=exp.label
	if (method=="Poisson")
	{
		colnames(results)=c("p", "lambda_S", "lambda_B", "k")
	}
	if (method=="NB")
	{
		colnames(results)=c("p", "mu_S", "phi_S", "mu_B", "phi_B", "k")
	}
	object=list(data=data, parameters=results, rep.vec=rep.vec, p.vec=p.vec, method=method)
	return(object)
}
