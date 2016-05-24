enrich.mix <-function(object, analysis="joint", differential=FALSE, diff.vec=NULL, cr=0.05, crdiff=0.05)
{
## INPUT
## object: the results of  mix if analysis="separate" or mix.joint if analysis="joint"
## differential: If TRUE, it will compute the posterior probability of differential binding of any two experiments or two conditions. 
## diff.vec: If differential = TRUE, diff.vec must be given to show which experiments are to be used in the comparison. At the moment, this is restricted to two conditions
## (e.g. two proteins at the same time point), so the value for diff.vec should be only 0, 1, 2, where, 0 indicates the experiments are not be used in the analysis, 
## 1 and 2 stand for conditions 1 and 2, respectively. diff.vec should be of the same length as the number of experiments in object.
## cr: the level of FDR for identifying the enriched regions.
## crdiff: the level of FDR for identifying the differentially bound regions.

	data=object$data
	datacount=as.matrix(data$count)
	Nexp=ncol(datacount)
	N=nrow(datacount)
	method=object$method
	results=object$parameters
	rep.vec=object$rep.vec
	p.vec=object$p.vec
	joint=ifelse(length(rep.vec)+length(p.vec)>0, 1, 0)
	exp.label=colnames(datacount)	
	if (is.null(exp.label)==1)
	{	
		exp.label=paste("Experiment", 1:Nexp, sep="")
	}
	ppx1=matrix(0, N, Nexp)
	colnames(ppx1)=exp.label
	ppx0=matrix(0, N, Nexp)
	colnames(ppx0)=exp.label
	X=matrix(0, N, Nexp)
	colnames(X)=exp.label
	IPEv=rep(0, Nexp)
	names(IPEv)=exp.label
	diffdata=NULL	
	diffprob1=NULL
	diffprob0=NULL
	diffX=NULL
	diffX1=NULL
	diffX2=NULL
	if (analysis=="joint"&joint==0)
		stop("You are using the joint method but you provide the separate modelling results. Please set analysis='separate'", call.=FALSE)
	if (analysis=="separate"&joint==1)
		stop("You are using the individual models for the joint modelling results. Please set analysis='joint'", call.=FALSE)
	if (analysis=="joint"&all(summary(factor(rep.vec))==1))
	{
		warning("No replicates used", call.=FALSE)
		analysis="separate"
	}		
	if (differential==TRUE&length(diff.vec)==0)
		stop("You must give the diff.vec vector, for detecting differentially bound regions", call.=FALSE)
	if (differential==TRUE&any(diff.vec>2)&any(diff.vec<0))
		stop("Only 0, 1, 2 can used in diff.vec: 0 means not in differential analysis, 1 means condition 1 and 2 means condition 2", call.=FALSE)
	for (i in 1:Nexp)
	{
		temp=pprob(datacount[,i], results[i,], method=method)
		ppx1[,i]=temp$px1
		ppx0[,i]=temp$px0
		temp1=FDR(temp$px0, cr=cr)
		X[,i]=temp1$X 
		IPEv[i]=IPE(results[i,], method=method)
	}
	if (analysis=="joint")
	{			
		if (length(rep.vec)>0)
		{
			rlist=factor(rep.vec)
			rlevel=levels(rlist)
			rlevelvalue=summary(rlist)
			srlevel=subset(rlevel, rlevelvalue>1&rlevel>0)
			srlevelvalue=subset(rlevelvalue, rlevelvalue>1&rlevel>0)
			srindex=matrix(0, Nexp, length(srlevelvalue))
			for (j in 1:length(srlevel))
			{
				srindex[,j]=ifelse(rlist==srlevel[j], 1, 0)
				data1=NULL
				para=NULL
				for (i in 1:Nexp)
				{
					if (srindex[i,j]==1)
					{
						data1=cbind(data1, datacount[,i])
						para=c(para, results[i,])
					}
				}
				pptemp=pprob_joint(data1, para, method=method)
				temp1=FDR(pptemp$px0, cr=cr)
				for (i in 1:Nexp)
				{
					if (srindex[i,j]==1)
					{
						ppx1[,i]=pptemp$px1
						ppx0[,i]=pptemp$px0
						X[,i]=temp1$X 
					}
				}
			}	
		}
	}
	if (differential==TRUE&all(diff.vec<=2&diff.vec>=0))
	{
		csindex=rep(0, Nexp)
		diffexp1=which(diff.vec==1)
		diffexp2=which(diff.vec==2)
		diffexp1.label=exp.label[diffexp1]
		diffexp2.label=exp.label[diffexp2]
		diffexp1.name="condition 1"
		diffexp2.name="condition 2"
		if (analysis=="joint")
		{
			csindex=rowSums(srindex)
		}
		diffindex=ifelse(diff.vec==1, -1, 0)+ifelse(diff.vec==2, 1, 0)
		diff1px1=rep(1, N)
		diff1px0=rep(1, N)
		diff2px1=rep(1, N)
		diff2px0=rep(1, N)
		if (sum(diffindex*csindex<0)>1)
		{
			data1=NULL
			para=NULL
			for (i in 1:Nexp)
			{
				if(diffindex[i]*csindex[i]==-1)
				{
					data1=cbind(data1, datacount[,i])
					para=c(para, results[i,])
				}				
			}
			temp=pprob_joint(data1, para, method=method)
			diff1px1=diff1px1*temp$px1
			diff1px0=diff1px0*temp$px0
			for (i in 1:Nexp)
			{
				if(diffindex[i]*csindex[i]>=0&diffindex[i]==-1)
				{
					diff1px1=diff1px1*ppx1[,i]
					diff1px0=diff1px0*ppx0[,i]
				}
			}
		}else
		{
			for (i in 1:Nexp)
			{
				if (diffindex[i]==-1)
				{
					diff1px1=diff1px1*ppx1[,i]
					diff1px0=diff1px0*ppx0[,i]
				}
			}
		}
		if (sum(diffindex*csindex>0)>1)
		{
			data1=NULL
			para=NULL
			for (i in 1:Nexp)
			{
				if(diffindex[i]*csindex[i]==1)
				{
					data1=cbind(data1, datacount[,i])
					para=c(para, results[i,])
				}				
			}
			temp=pprob_joint(data1, para, method=method)
			diff2px1=temp$px1
			diff2px0=temp$px0
			for (i in 1:Nexp)
			{
				if(diffindex[i]*csindex[i]>=0&diffindex[i]==1)
				{
					diff2px1=diff2px1*ppx1[,i]
					diff2px0=diff2px0*ppx0[,i]
				}
			}
		}else
		{
			for (i in 1:Nexp)
			{
				if (diffindex[i]==1)
				{
					diff2px1=diff2px1*ppx1[,i]
					diff2px0=diff2px0*ppx0[,i]
				}
			}
		}
		diffprob1=diff1px1*diff2px0+diff1px0*diff2px1 
		diffprob0=1-diffprob1
		temp1<-FDR(diffprob0, cr=crdiff)
		diffX<-temp1$X
		diffX1<-ifelse(diffX==1&diff1px1>diff2px1, 1, 0)
		diffX2<-ifelse(diffX==1&diff1px1<diff2px1, 1, 0)
		diffX1=as.matrix(diffX1)
		colnames(diffX1)=diffexp1.name
		diffX2=as.matrix(diffX2)
		colnames(diffX2)=diffexp2.name
	}

	## enrich and differential bound regions
	nenrich=colSums(X)
	enrich<-list()
	diffenrich1<-NULL
	diffenrich2<-NULL
	cat ("***************Enrich regions**************", "\n")
	cat ("There are", Nexp, "experiments.", '\n')
	cat ("The number of enriched regions and IP efficiency estimates for these experiments are:", '\n')
	rlist=factor(rep.vec)
	rlevel=levels(rlist)
	rlevelvalue=summary(rlist)
	srlevel=subset(rlevel, rlevelvalue>1&rlevel>0)
	srlevelvalue=subset(rlevelvalue, rlevelvalue>1&rlevel>0)
	flist=factor(p.vec)
	level=levels(flist)
	levelvalue=summary(flist)
	splevel=subset(level, levelvalue>1)
	splevelvalue=subset(levelvalue, levelvalue>1)
	spindex=matrix(0, Nexp, length(splevelvalue))
	if (joint==1&length(splevel)>0)
	{	
		flist=factor(p.vec)
		level=levels(flist)
		levelvalue=summary(flist)
		splevel=subset(level, levelvalue>1)
		splevelvalue=subset(levelvalue, levelvalue>1)
		spindex=matrix(0, Nexp, length(splevelvalue))
		for (i in 1:length(splevel))
		{
			spindex[,i]=ifelse(flist==splevel[i], 1, 0)
			spcode=which(spindex[,i]==1)
			templabel=exp.label[spcode]
			cat("Note: in previous analysis", templabel, "are assumed have same p", '\n')
		}
	}
	if(joint==1&length(srlevel)>0)
	{
		for (i in 1:length(srlevel))
		{
			jointcode=which(srindex[,i]==1)
			templabel=exp.label[jointcode]
			tempnenrich=nenrich[jointcode][1]
			tempIPEv=IPEv[jointcode]
			if(length(jointcode)==2)
			{
				templabel=paste(templabel[1], "and", templabel[2], sep=" ")
				tempIPEv=paste(tempIPEv[1], "and", tempIPEv[2], sep=" ")
			}
			cat(templabel, " are replicates and they have", tempnenrich, "jointly identified enriched regions at FDR=", cr, "and their IP efficiencies are ", tempIPEv, "respectively", '\n')
			temp=subset(cbind(data$region, data$count[,jointcode], ppx1[,jointcode[1]]), X[,jointcode[1]]==1)
			colnames(temp)=c(colnames(data$region), exp.label[jointcode], "poster.prob")
			temp=temp[order(temp[,1], temp[,2]),]
			enrich[[i]]=temp	
		}
		sepcode=which(rowSums(srindex)==0)
		if (length(sepcode)>0)
		{
			for (j in 1:length(sepcode))
			{
				i=sepcode[j]
				cat (exp.label[i], "has", nenrich[i], "enriched regions at FDR=", cr, "and its IP efficiency is", IPEv[i], '\n')
				temp=subset(cbind(data$region, data$count[,i], ppx1[,i]), X[,i]==1)
				colnames(temp)=c(colnames(data$region), exp.label[i], "poster.prob")
				temp=temp[order(temp[,1], temp[,2]),]
				enrich[[j+length(srlevel)]]=temp	
			}
		}
	}else
	{
		for (i in 1:Nexp)
		{
			cat (exp.label[i], "has", nenrich[i], "enriched regions at FDR=", cr, "and its IP efficiency is", IPEv[i], '\n')
			temp=subset(cbind(data$region, data$count[,i], ppx1[,i]), X[,i]==1)
			colnames(temp)=c(colnames(data$region), exp.label[i], "poster.prob")
			temp=temp[order(temp[,1], temp[,2]),]
			enrich[[i]]=temp	
		}
	}
	if (length(diffX1)>0)
	{
		cat ("***************Differential bound regions**************", "\n")
		cat ("Condition 1 uses data of", diffexp1.label, '\n')
		cat ("Condition 2 uses data of", diffexp2.label, '\n')
		cat ("The number of regions bound only by condition 1 at FDR=", crdiff, " is:", sum(diffX1), '\n')
		cat ("The number of regions bound only by condition 2 at FDR=", crdiff, " is:", sum(diffX2), '\n')
		diffenrich1<-subset(cbind(data$region, data$count[,diffexp1], diffprob1), diffX1==1)
		colnames(diffenrich1)=c(colnames(data$region), diffexp1.label, "poster.prob.diff")
		diffenrich2<-subset(cbind(data$region, data$count[,diffexp2], diffprob1),  diffX2==1)
		colnames(diffenrich2)=c(colnames(data$region), diffexp2.label, "poster.prob.diff")
	}

##    Output list
##    enrich: the list of enrich regions. 
##    diffenrich1 and diffenrich2: the differential bound regions of the first and second conditions, respectively.
##    ppx1: n x p matrix of posterior probabilities of each region and experiment. ppx0=1-ppx1
##    X: n x p matrix of index of enrichment for each experiment (at the FDR cutoff), X=1 : enriched region, X=0 : background region
##    IPE: p-dimensional vector of IP efficiency values for the p experiments.
##    diffprob1: n-dimensional vector of posterior probabilities P(X_1 \ne X_2) for the two conditions under study; diffprob0=1-diffprob1
##    diffX1: n-dimensional vector of index of differential bound regions for condition 1, diffX2 is the index of differential bound regions for condition 2.
 
	object<-list(enrich=enrich, diffenrich1=diffenrich1, diffenrich2=diffenrich2, ppx1=ppx1, ppx0=ppx0, X=X, IPE=IPEv, diffprob1=diffprob1, diffprob0=diffprob0, diffX1=diffX1, diffX2=diffX2)
	return(object)
}
