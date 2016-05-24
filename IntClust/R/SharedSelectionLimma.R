SharedSelectionLimma<-function(DataLimma=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	
	which=list()	
	table=c()
	
	if(is.null(names)){
		for(j in 1:length(DataLimma)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	nmethods=length(DataLimma)
	
	temp1g=c()
	comps=c()
	
	pvalsg=c()
	for (i in 1:nmethods){			
		
		temp1g=c(temp1g,length(DataLimma[[i]]$Genes$TopDE$Genes))
		comps=c(comps,length(DataLimma[[i]]$Compounds$LeadCpds))
		
		
		names(temp1g)[i]=names[i]
		names(comps)[i]=paste("Ncomps",names[i],i,sep=" ")
		
		
		
		
		if (i==1){
			if(!(is.na(DataLimma[[i]])[1])){
				sharedcomps=DataLimma[[i]]$Compounds$LeadCpds
				sharedgenes=DataLimma[[i]]$Genes$TopDE$Genes
				
				
				pvalsg=c(pvalsg,DataLimma[[i]]$Genes$TopDE$adj.P.Val)
				
				
				nsharedcomps=length(DataLimma[[i]]$Compounds$LeadCpds)
				nsharedgenes=length(DataLimma[[i]]$Genes$TopDE$Genes)
				
				names(nsharedgenes)="nshared"
				
				names(nsharedcomps)="nsharedcomps"
			}
			
		}
		else{			
			sharedcomps=intersect(sharedcomps,DataLimma[[i]]$Compounds$LeadCpds)
			sharedgenes=intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$Genes)
			
			
			nsharedcomps=length(intersect(sharedcomps,DataLimma[[i]]$Compounds$LeadCpds))
			nsharedgenes=length(intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$Genes))
			
			names(nsharedgenes)="nshared"
			
			names(nsharedcomps)="nsharedcomps"
		}
		
	}	
	pvalsgenes=list()
	meanpvalsgenes=c()
	
	if(nsharedgenes != 0){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalsg=c(pvalsg,DataLimma[[c]]$Genes$TopDE$adj.P.Val[DataLimma[[c]]$Genes$TopDE$Genes==g])	
				}	
			}
			
			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("Method",c,sep=" ")
		}	
		
		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()			
			for(c in 1:nmethods){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}			
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes	
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{pvalsgenes=0}
	
	
	temp=rbind(temp1g,nsharedgenes,comps,nsharedcomps)	
	
	table=cbind(table,temp)
	
	which[[1]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes)
	#names(which)[1]=paste("Cluster",i,sep=" ")
	
	
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)	
}
