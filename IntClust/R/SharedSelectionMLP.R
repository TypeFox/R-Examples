SharedSelectionMLP<-function(DataMLP=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	
	which=list()	
	table=c()
	
	
	nmethods=length(DataMLP)
	
	if(is.null(names)){
		for(j in 1:length(DataMLP)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	temp1g=c()
	temp1p=c()
	comps=c()
	
	pvalsg=c()
	pvalsp=c()	
	for (i in 1:nmethods){			
		
		temp1g=c(temp1g,length(DataMLP[[i]]$Genes$TopDE$Genes))
		temp1p=c(temp1p,length(DataMLP[[i]][[3]]$geneSetDescription))
		comps=c(comps,length(DataMLP[[i]]$Compounds$LeadCpds))
		
		
		names(temp1g)[i]=names[i]
		names(temp1p)[i]=names[i]
		names(comps)[i]=paste("Ncomps",names[i],i,sep=" ")
		
		
		if (i==1){
			if(!(is.na(DataMLP[[i]])[1]) | !(is.na(DataMLP[[i]])[1])){
				sharedcomps=DataMLP[[i]]$Compounds$LeadCpds
				sharedgenes=DataMLP[[i]]$Genes$TopDE$Genes
				sharedpaths=DataMLP[[i]][[3]]$geneSetDescription
				
				pvalsg=c(pvalsg,DataMLP[[i]]$Genes$TopDE$adj.P.Val)
				pvalsp=c(pvalsp,DataMLP[[i]]$mean_geneSetPValue)
				
				nsharedcomps=length(DataMLP[[i]]$Compounds$LeadCpds)
				nsharedgenes=length(DataMLP[[i]]$Genes$TopDE$Genes)
				nsharedpaths=length(DataMLP[[i]][[3]]$geneSetDescription)
				names(nsharedgenes)="nshared"
				names(nsharedpaths)="nshared"
				names(nsharedcomps)="nsharedcomps"
			}
			
		}
		else{			
			sharedcomps=intersect(sharedcomps,DataMLP[[i]]$Compounds$LeadCpds)
			sharedgenes=intersect(sharedgenes,DataMLP[[i]]$Genes$TopDE$Genes)
			sharedpaths=intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription)
			
			nsharedcomps=length(intersect(sharedcomps,DataMLP[[i]]$Compounds$LeadCpds))
			nsharedgenes=length(intersect(sharedgenes,DataMLP[[i]]$Genes$TopDE$Genes))
			nsharedpaths=length(intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription))
			names(nsharedgenes)="nshared"
			names(nsharedpaths)="nshared"
			names(nsharedcomps)="nsharedcomps"
		}
		
	}	
	pvalsgenes=list()
	meanpvalsgenes=c()
	meanpvalspaths=c()
	pvalspaths=list()
	
	if(nsharedgenes != 0){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalsg=c(pvalsg,DataMLP[[c]]$Genes$TopDE$adj.P.Val[DataMLP[[c]]$Genes$TopDE$Genes==g])	
				}	
			}
			
			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("Method",c,sep=" ")
		}	
		
		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()			
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}			
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes	
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{pvalsgenes=0}
	
	if(nsharedpaths!=0){
		for(c in 1:nmethods){
			pvalsp=c()
			if(!(is.na(DataMLP[[c]])[1])){
				for(p in sharedpaths){
					pvalsp=c(pvalsp,DataMLP[[c]][[3]][DataMLP[[c]][[3]]$geneSetDescription==p,5])
				}
			}
			
			pvalspaths[[c]]=pvalsp
			
			
			names(pvalspaths)[c]=paste("Method",c,sep=" ")
		}
		
		
		for(p1 in 1:length(sharedpaths)){
			pvalstemp1=c()
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])
					
				}
				
			}			
			
			meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
			
		}
		pvalspaths[[nmethods+1]]=meanpvalspaths	
		names(pvalspaths)[nmethods+1]="Mean pvals paths"
	}
	else{pvalpaths=0}
	
	temp=rbind(cbind(temp1g,temp1p),cbind(nsharedgenes,nsharedpaths),cbind(comps,comps),cbind(nsharedcomps,nsharedcomps))	
	
	table=cbind(table,temp)
	
	which[[1]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes,sharedpaths=sharedpaths,pvalspaths=pvalspaths)
	#names(which)[1]=paste("Cluster",i,sep=" ")
	
	
	for (i in 1:length(seq(1,dim(table)[2],2))){
		number=seq(1,dim(table)[2],2)[i]
		colnames(table)[number]=c("G.Cluster")
		colnames(table)[number+1]=c("P.Cluster")	
		
	}
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)
	
	
}
