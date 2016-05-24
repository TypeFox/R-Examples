Geneset.intersect<-function(PathwaysOutput,Selection=FALSE,sign=0.05,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	
	if(Selection==TRUE){
		if(length(PathwaysOutput$'Iteration 1')==1 & is.null(names)){
			names="Selection"
		}
		Intersect=Geneset.intersectSelection(PathwaysOutput,sign,names,seperatetables,separatepvals)	
	}
	
	else{
	
	if(is.null(names)){
		for(j in 1:length(PathwaysOutput$"Iteration 1")){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	
	#put all of same method together:preparation of lists
	subsets=list()
	nmethods=length(PathwaysOutput$"Iteration 1") 
	for(i in 1:nmethods){
		subsets[[i]]=list()
		
	}
	names(subsets)=names
	
	#put all of same method together: go through PathwaysOutput
	for(j in 1:length(PathwaysOutput)){
		name1=names(PathwaysOutput)[j]
		for(k in 1:nmethods){
			name2=names[k]
			subsets[[name2]][[name1]]=PathwaysOutput[[name1]][[name2]]
			
		}
		
	}	
	
	#for every subset (= every method) take intersection over the interations per cluster
	Intersect=list()
	
	for(i in 1:length(subsets)){
		Method=subsets[[i]]
		Clusters=list()
		nclus=length(Method[[1]])
		for(j in 1:length(Method)){
			name3=paste("Iteration",j,sep=" ")
			for(k in 1:nclus){
				name4=paste("Cluster",k,sep=" ")
				if(!(is.na(Method[[name3]][[name4]])[1])){
					Clusters[[name4]][[name3]]=Method[[name3]][[name4]]
				}
				else{
					Clusters[[name4]][[name3]]=NA
				}
			}
			
		}
		
		IntersectM=list()
		
		for(a in 1:length(Clusters)){ #per cluster
			if(!(is.na(Clusters[[a]])[1])){
				result.out=list()
				result.name = c()
				for(b in 1:length(Clusters[[a]])){#per iteration
					if(b==1){
						Compounds=Clusters[[a]][[1]]$Compounds
						Genes=Clusters[[a]][[1]]$Genes	
						Names=data.frame("description"=Clusters[[a]][[1]]$Pathways$AllPaths$geneSetDescription,"genesetcode"=rownames(Clusters[[a]][[1]]$Pathways$AllPaths))
						Names$description=as.character(Names$description)	
						Names$genesetcode=as.character(Names$genesetcode)	
					}				
					cut = Clusters[[a]][[b]]$Pathways$AllPaths[  Clusters[[a]][[b]]$Pathways$AllPaths[,2]<=sign,]
					colnames(cut)[4] = paste("pvalues.",b,sep="")
					colnames(cut)[2] = paste("testedgenesetsize.",b,sep="")
					colnames(cut)[3] = paste("genesetstatistic.",b,sep="")
					cut=cut[,c(1,5,2,3,4)]
					result.out[[b]] = cut
					result.name = c(result.name,paste("genesettable",b,sep=""))

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
				
				
				
				if(separatepvals==FALSE){
					result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,(ncol(genesets.table.intersect)-2):ncol(genesets.table.intersect))]
					class(result.out$genesets.table.intersect)=c("MLP","data.frame")
					attr(result.out$genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
				}
				
				
				if(seperatetables==FALSE){
					result.out=result.out$genesets.table.intersect
					class(result.out)=c("MLP","data.frame")
					attr(result.out)=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
				}
				
				newresult=list(Compounds=Compounds,Genes=Genes,Pathways=result.out)
				
				
				IntersectM[[a]]=newresult	
				names(IntersectM)[a]=names(Clusters)[[a]]
			}
			else{
				IntersectM[[a]]=NA	
				names(IntersectM)[a]=names(Clusters)[[a]]
			}
		}
		
		Intersect[[i]]=IntersectM
		
		}
	}
	names(Intersect)=names
	return(Intersect)
}
