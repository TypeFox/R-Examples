SharedGenesPathsFeat<-function(DataLimma=NULL,DataMLP=NULL,DataFeat=NULL,names=NULL,Selection=FALSE){  #Input=result of DiffGenes and Geneset.intersect
	#Include sharedLimma and SharedMLP inside the function
	if(Selection==TRUE){
		ResultShared=SharedSelection(DataLimma,DataMLP,DataFeat,names)
	}	
	
	else if(is.null(DataLimma) & is.null(DataMLP) & is.null(DataFeat)){	
		stop("At least one Data set should be specified")
	}
	
	else{
	
	List=list(DataLimma,DataMLP,DataFeat)
	AvailableData=sapply(seq(length(List)),function(i) if(!(is.null(List[[i]]))) return(i))
	AvailableData=unlist(AvailableData)
	len=c()
	for(i in AvailableData){
	 	len=c(len,length(List[[i]]))
	}
	if(length(unique(len))!=1){
	 stop("Unequal number of methods for limma and MLP")
	}
	else{
	 DataSets=lapply(AvailableData,function(i)  return(List[[i]]))
	 nmethods=length(DataSets[[1]])
	 nclusters=length(DataSets[[1]][[1]])
	}
	
	if(is.null(names)){
		for(j in 1:length(DataSets[[1]])){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	which=list()	
	table=c()
		
			
	for (i in 1:nclusters){
		
		name=paste("Cluster",i,sep=" ")
		
		comps=c()
		
		temp1g=c()
		temp1p=c()
		temp1f=list()

		
		for(j in 1:nmethods){

			if(!(is.na(DataSets[[1]][[j]][[i]])[1])){
				comps=c(comps,length(DataSets[[1]][[j]][[i]]$Compounds$LeadCpds))
				names(comps)[j]=paste("Ncomps", names[j],sep=" ")
			}
			else{
				comps=c(comps,"-")
			}
		
			if(!(is.null(DataLimma))){
				if(!(is.na(DataLimma[[j]][[i]])[1])){
					temp1g=c(temp1g,length(DataLimma[[j]][[i]]$Genes$TopDE$ID))
				}
				else{
					temp1g=c(temp1g,"-")
				}	
				names(temp1g)[j]=names[j]
			}
			else{
				temp1g=NULL
			}
		
			if(!(is.null(DataMLP))){
				if(!(is.na(DataMLP[[j]][[i]])[1])){
					temp1p=c(temp1p,length(DataMLP[[j]][[i]][[3]]$geneSetDescription))
				}
				else{
					temp1p=c(temp1g,"-")
				}	
				names(temp1p)[j]=names[j]
			}
			else{
				temp1p=NULL
			}
			
			if(!(is.null(DataFeat))){
				temp=c()
				for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){			
					if(!(is.na(DataFeat[[j]][[i]])[1])){
						temp=c(temp,length(DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names))
					}
					else{
						temp=c(temp,"-")
					}
					
					names(temp)[f]=names(DataFeat[[j]][[i]]$Characteristics)[f]
				}	
				temp1f[[j]]=temp
				names(temp1f)[j]=names[j]	
			}
			else{
				temp1f=NULL
			}
		}	

		j=1		
		Continue=TRUE
		while (Continue==TRUE){
			cont=c()
			for(d in 1:length(DataSets)){
				cont=c(cont,!(is.na(DataSets[[d]][[j]][[i]])[1]))
			}	

			if(any(cont)){
				
				sharedcomps=DataSets[[1]][[j]][[i]]$Compounds$LeadCpds
				nsharedcomps=length(sharedcomps)
				names(nsharedcomps)="nsharedcomps"
			
			
				if(!(is.null(DataLimma))){
					sharedgenes=DataLimma[[j]][[i]]$Genes$TopDE$ID
					nsharedgenes=length(sharedgenes)
					names(nsharedgenes)="Nshared"
					#pvalsg=DataLimma[[j]][[i]]$Genes$TopDE$adj.P.Val
				}
				else{
					sharedgenes=NULL
					nsharedgenes=NULL
				}
					

				if(!(is.null(DataMLP))){
					sharedpaths=DataMLP[[j]][[i]][[3]]$geneSetDescription
					nsharedpaths=length(sharedpaths)
					names(nsharedpaths)="Nshared"
					#pvalsp=DataMLP[[j]][[i]][[3]]$mean_geneSetPValue
					
				}
				else{
					sharedpaths=NULL
					nsharedpaths=NULL
				}
				if(!(is.null(DataFeat))){
					sharedfeat=list()
					nsharedfeat=list()
					#pvalsf=list()
					for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){	
						sharedfeat[[f]]=DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names	
						names(sharedfeat)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
						
						nsharedfeat[[f]]=length(sharedfeat[[f]])
						names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
						
						#pvalsf[[f]]=DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$adj.P.Val	
						#names(pvalsf)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
						
					}
								
				}
				else{
					sharedfeat=NULL
					nsharedfeat=NULL
				}
			
				Continue=FALSE
			}
			j=j+1
		}

		if(nmethods>=2){
			for (j in 2:nmethods){
				cont=c()
				for(d in 1:length(DataSets)){
					cont=c(cont,!(is.na(DataSets[[d]][[j]][[i]])[1]))
				}	
				if(any(cont)){
					
					sharedcomps=intersect(sharedcomps,DataSets[[1]][[j]][[i]]$Compounds$LeadCpds)
					nsharedcomps=length(sharedcomps)
					names(nsharedcomps)="Nsharedcomps"
					
					
					if(!(is.null(DataLimma))){
						sharedgenes=intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$TopDE$ID)
						nsharedgenes=length(sharedgenes)
						names(nsharedgenes)="Nshared"
						
					}
					if(!(is.null(DataMLP))){
						sharedpaths=intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$geneSetDescription)
						nsharedpaths=length(sharedpaths)
						names(nsharedpaths)="Nshared"
						
					}
					if(!(is.null(DataFeat))){

						for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){	
							sharedfeat[[f]]=intersect(sharedfeat[[f]],DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names)	
							names(sharedfeat)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
							
							nsharedfeat[[f]]=length(sharedfeat[[f]])
							names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
							
						}
						
					}
			
				}
			}			
		}

		pvalsgenes=list()
		meanpvalsgenes=c()
		
		pvalspaths=list()
		meanpvalspaths=c()
		
		pvalsfeat=list()
		meanpvalsfeat=c()
		
			
		if(!(is.null(sharedgenes))&length(sharedgenes)!=0){
			for(c in 1:nmethods){
				pvalsg=c()
				for(g in sharedgenes){
					if(!(is.na(DataLimma[[c]][[i]])[1])){
						pvalsg=c(pvalsg,DataLimma[[c]][[i]]$Genes$TopDE$adj.P.Val[DataLimma[[c]][[i]]$Genes$TopDE$ID==g])	
					}	
				}
				
				pvalsgenes[[c]]=pvalsg
				names(pvalsgenes)[c]=paste("P.Val.",names[c],sep="")
			}	
				
			for(g1 in 1:length(sharedgenes)){
				pvalstemp=c()			
				for(c in 1:nmethods){
					if(!(is.na(DataLimma[[c]][[i]])[1])){
						pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
					}
				}			
				meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
			}
			pvalsgenes[[nmethods+1]]=meanpvalsgenes	
			names(pvalsgenes)[nmethods+1]="Mean pvals genes"
		}
		else{pvalsgenes=NULL}
					
		if(!(is.null(sharedpaths))&length(sharedpaths)!=0){
			for(c in 1:nmethods){
				pvalsp=c()
				if(!(is.na(DataMLP[[c]][[i]])[1])){
					for(p in sharedpaths){
						pvalsp=c(pvalsp,DataMLP[[c]][[i]][[3]][DataMLP[[c]][[i]][[3]]$geneSetDescriptions==p,5][1])
					}
				}
					
				pvalspaths[[c]]=pvalsp
				names(pvalspaths)[c]=paste("P.Val.",names[c],sep="")
			}
				
				
			for(p1 in 1:length(sharedpaths)){
				pvalstemp1=c()
				for(c in 1:nmethods){
					if(!(is.na(DataMLP[[c]][[i]])[1])){
						pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])
						
					}
				
				}			
					
				meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
			}
			pvalspaths[[nmethods+1]]=meanpvalspaths	
			names(pvalspaths)[nmethods+1]="Mean pvals paths"
		}
		else{pvalpaths=NULL}
		
		
		if(!(is.null(sharedfeat))){
			for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){

				if(length(sharedfeat[[f]])!=0){
					pvalschar=list()
					for(c in 1:nmethods){
						pvalsf=c()
						if(!(is.na(DataFeat[[c]][[i]])[1])){
							for(s in sharedfeat[[f]]){
								pvalsf=c(pvalsf,DataFeat[[c]][[i]]$Characteristics[[f]]$TopFeat$adj.P.Val[DataFeat[[c]][[i]]$Characteristics[[f]]$TopFeat$Names==s])
							}
						}
						pvalschar[[c]]=pvalsf
						names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
					}
				}
				else{
					pvalschar=list()
					pvalschar[1:nmethods]=0
					for(c in 1:nmethods){
						names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
					}
				}
				pvalsfeat[[f]]=pvalschar
				names(pvalsfeat)[f]=names(DataFeat[[j]][[i]]$Characteristics)[f]

			}

			for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
				meanpvalsfeat=c()
	
				if((length(sharedfeat[[f]])!=0)){				
					for(f1 in 1:length(sharedfeat[[f]])) {
	
						pvalstemp=c()			
						for(c in 1:nmethods){
							if(!(is.na(DataFeat[[c]][[i]])[1])){
								pvalstemp=c(pvalstemp,pvalsfeat[[f]][[c]][[f1]])
							}
						}			
						meanpvalsfeat=c(meanpvalsfeat,mean(pvalstemp))			
					}
				}
				else{
					meanpvalsfeat=0
				}
				
				pvalsfeat[[f]][[nmethods+1]]=meanpvalsfeat
				names(pvalsfeat[[f]])[nmethods+1]="Mean pvals feat"
			}
			lenchar=length(DataFeat[[j]][[i]]$Characteristics)
			
		}
		else{
			pvalsfeat=NULL
			lenchar=0
		}
		
			
		if(!(is.null(temp1f))){
			temp1f=do.call(rbind.data.frame, temp1f)	
			colnames(temp1f)=names(DataFeat[[j]][[i]]$Characteristics)
			temp1f=as.matrix(temp1f)
			nsharedfeat=do.call(cbind.data.frame, nsharedfeat)
			nsharedfeat=as.matrix(nsharedfeat)
		}
		part1=cbind(cbind(temp1g,temp1p),temp1f)
		part1=as.matrix(part1)
		if(is.null(nsharedgenes) & is.null(nsharedpaths) &is.null(nsharedfeat)){
			if(!(is.null(temp1g))){
				nsharedgenes=0
			}
			if(!(is.null(temp1p))){
				nsharedpaths=0
			}
			if(!(is.null(temp1f))){
				nsharedfeat=rep(0,length(temp1f))
			}
		}
		part2=cbind(cbind(nsharedgenes,nsharedpaths),nsharedfeat)
		part2=as.matrix(part2)
		colnames(part1)=NULL
		colnames(part2)=NULL
		rownames(part2)="NShared"
		part3=c()
		for(r in 1:length(comps)){
			part3=rbind(part3,rep(comps[r],dim(part1)[2]))
		}
		colnames(part3)=NULL
		rownames(part3)=names(comps)
		part4=rep(nsharedcomps,dim(part1)[2])
		names(part4)=NULL
		temp=rbind(part1,part2,part3,part4)	
		rownames(temp)[nrow(temp)]="Nsharedcomps"

		
		table=cbind(table,temp)
		colnames(table)=seq(1,ncol(table))
		
		if(!(is.null(pvalsgenes))){
			SharedGenes=cbind(sharedgenes,do.call(cbind.data.frame, pvalsgenes))
		}
		else{
			SharedGenes=NULL
		}
		if(!(is.null(pvalspaths))){
			SharedPaths=cbind(sharedpaths,do.call(cbind.data.frame, pvalspaths))
		}
		else{
			SharedPaths=NULL
		}
		if(!(is.null(pvalsfeat))){
			SharedFeat=list()
			for(f in 1:lenchar){
				if(length(sharedfeat[[f]])==0){
					SharedFeat[[f]]=NULL
				}
				else{
					SharedFeat[[f]]=cbind(sharedfeat[[f]],do.call(cbind.data.frame, pvalsfeat[[f]]))
					names(SharedFeat)[f]=names(pvalsfeat[[1]])[f]
				}
				
			}
		}
		else{
			SharedFeat=NULL
		}
		which[[i]]=list(SharedComps=sharedcomps,SharedGenes=SharedGenes,SharedPaths=SharedPaths,SharedFeat=SharedFeat)
		names(which)[i]=paste("Cluster",i,sep=" ")
		
	}
	#Sep for all situations?
	if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+2))){
			number=seq(1,dim(table)[2],2+lenchar)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
			colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")	
			for(u in seq(1:lenchar)){
				colnames(table)[number+1+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}
		
		}
	}	
	
	else if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),is.null(DataFeat))){
		for (i in 1:length(seq(1,dim(table)[2],2))){
			number=seq(1,dim(table)[2],2)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
			colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")	
			
		}
	}	
	
	else if(all(!(is.null(DataLimma)),is.null(DataMLP),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
			number=seq(1,dim(table)[2],1+lenchar)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
			for(u in seq(1:lenchar)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}
			
		}
	}
	
	else if(all((is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
			number=seq(1,dim(table)[2],lenchar+1)[i]	
			colnames(table)[number]=paste("P.Cluster",i,sep=" ")	
			for(u in seq(1:lenchar)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}
			
		}
	}	
	
	else if(all(!(is.null(DataLimma)),(is.null(DataMLP)),(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],1))){
			colnames(table)[i]=paste("G.Cluster",i,sep=" ")				
		}
	}	
	
	else if(all((is.null(DataLimma)),!(is.null(DataMLP)),(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],1))){
			colnames(table)[i]=paste("P.Cluster",i,sep=" ")				
		}
	}	
	
	else if(all((is.null(DataLimma)),(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar))){
			number=seq(1,dim(table)[2],2)[i]
			for(u in c(0,1)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u+1],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}
			
		}
	}	
	
	ResultShared=list(Table=table,Which=which)	
	
	}
	return(ResultShared)
	
}
