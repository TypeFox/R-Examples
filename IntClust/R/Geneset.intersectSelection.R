Geneset.intersectSelection<-function(list.output,sign,names=NULL,seperatetables=FALSE,separatestats=FALSE,geneSetSource = "GOBP"){
	if(is.null(names)){
		for(j in 1:length(list.output$"Iteration 1")){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	#put all of same method together:preparation of lists
	subsets=list()
	nmethods=length(list.output$"Iteration 1") 
	for(i in 1:nmethods){
		subsets[[i]]=list()		
	}
	names(subsets)=names

	#put all of same method together: go through list.output
	for(j in 1:length(list.output)){
		name1=names(list.output)[j]
		for(k in 1:nmethods){
			name2=k
			subsets[[name2]][[name1]]=list.output[[name1]][[name2]]
			
		}
		
	}	
	
	#for every subset (= every method) take intersection over the interations per cluster
	Intersect=list()
	
	for(i in 1:length(subsets)){
		Method=subsets[[i]]
		Clusters=list()
		nclus=1
		for(j in 1:length(Method)){
			name3=paste("Iteration",j,sep=" ")
			Clusters[[name3]]=Method[[name3]]			
		}
		
		IntersectM=list()
		
		result.out=list()
		result.name = c()
		for(a in 1:length(Clusters)){ #per cluster
			if(a==1){
				Compounds=Clusters[[a]]$Compounds
				Genes=Clusters[[a]]$Genes
				Names=data.frame("description"=Clusters[[a]]$Pathways$AllPaths$geneSetDescription,"genesetcode"=rownames(Clusters[[a]]$Pathways$AllPaths))
				Names$description=as.character(Names$description)	
				Names$genesetcode=as.character(Names$genesetcode)	
			}				
			cut = Clusters[[a]]$Pathways$AllPaths[Clusters[[a]]$Pathways$AllPaths$geneSetPValue<=sign,]
			colnames(cut)[4] = paste("pvalues.",a,sep="")
			colnames(cut)[2] = paste("testedgenesetsize.",a,sep="")
			colnames(cut)[3] = paste("genesetstatistic.",a,sep="")
			cut=cut[,c(1,5,2,3,4)]
			result.out[[a]] = cut
			result.name = c(result.name,paste("genesettable",a,sep=""))
		}
			
		names(result.out) = result.name
		
		genesets.table.intersect = plyr::join_all(result.out,by=c("totalGeneSetSize","geneSetDescription"),type="inner")
		genesets.table.intersect$mean_testedGeneSetSize=round(apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='testedgenesetsize')],1,mean),1)
		genesets.table.intersect$mean_geneSetStatistic=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='genesetstatistic')],1,mean)
		genesets.table.intersect$mean_geneSetPValue=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='pvalues')],1,mean)
		
		rownames(genesets.table.intersect)=as.character(Names[which(genesets.table.intersect$geneSetDescription%in%Names[,1]),2])
		
		class(genesets.table.intersect)=c("MLP","data.frame")
		attr(genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource

		
		result.out$genesets.table.intersect = genesets.table.intersect
		
		if(separatestats==FALSE){
			result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,(ncol(genesets.table.intersect)-2):ncol(genesets.table.intersect))]
			class(result.out$genesets.table.intersect)=c("MLP","data.frame")
			attr(result.out$genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		}
		
		
		if(seperatetables==FALSE){
			result.out=result.out$genesets.table.intersect
			class(result.out)=c("MLP","data.frame")
			attr(result.out,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		}
		
		
		newresult=list(Compounds=Compounds,Genes=Genes,Pathways=result.out)
		
		
		#IntersectM[[a]]=
		#names(IntersectM)[a]=names(Clusters)[[a]]
		Intersect[[i]]=newresult	
		
	}
	names(Intersect)=names
	return(Intersect)
}
