Pathways<-function(List,Selection=NULL,GeneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),GeneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	
	
	if(!(is.null(Selection))){
		ResultMLP=PathwaysSelection(List,Selection,GeneExpr,nrclusters,method,GeneInfo,geneSetSource,topP,topG,GENESET,sign,fusionsLog,WeightClust,names)
			
	}
	else if(class(List)=="ChosenClusters"){
		ResultMLP=list()
		for(i in 1:length(List)){	
			Selection=List[[i]]$Compounds$LeadCpds
			L=List[i]
			ResultMLP[[i]]=PathwaysSelection(List=L,Selection,GeneExpr,nrclusters,method,GeneInfo,geneSetSource,topP,topG,GENESET,sign,fusionsLog,WeightClust,names)	
			names(ResultMLP)=paste("Choice",i,sep=' ')
		}	
	}
	else{
		
		#Check for gene expression data: if not, reordering to ListNew not necessary
		DataPrepared<-plyr::try_default(PreparePathway(List[[1]],GeneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(DataPrepared)){
		
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(WeightClust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		List=ListNew	
		
		maxclus=0
		DataPrepared=list()
		for (k in 1:dim(MatrixClusters)[1]){
			message(k)
			clusters=MatrixClusters[k,]
			
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			
			check<-plyr::try_default(PreparePathway(List[[k]],GeneExpr,topG,sign),NULL,quiet=TRUE)
			if(is.null(check)){
				Temp=List[[k]]
				
				for(i in unique(clusters)){
					Compounds=list()
					Compounds$LeadCpds=names(clusters)[which(clusters==i)] 
					Compounds$OrderedCpds=stats::as.hclust(List[[k]]$Clust$Clust)$labels[stats::as.hclust(List[[k]]$Clust$Clust)$order]
					
					Temp[[i+1]]=list(Compounds=Compounds)
					names(Temp)[i+1]=paste("Cluster",i,sep=" ")
					
				}
				DataPrepared[[k]]<-PreparePathway(Temp,GeneExpr,topG,sign)
				
			}	
		}
		}
		
		else{
			for(k in 1:length(List)){
				DataPrepared[[k]]<-plyr::try_default(PreparePathway(List[[k]],GeneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(DataPrepared[[k]])){
					Temp=List[[k]]
					
					for(i in unique(clusters)){
						Compounds=list()
						Compounds$LeadCpds=List[[k]]$Compounds$LeadCpds
						Compounds$OrderedCpds=List[[k]]$Compounds$OrderedCpds
						
						Temp[[i+1]]=list(Compounds=Compounds)
						names(Temp)[i+1]=paste("Cluster",i,sep=" ")
						
					}
					DataPrepared[[k]]<-PreparePathway(Temp,GeneExpr,topG,sign)
					
				}	
			}
			
			
		}
		
		

		method.test = function(sign.method,path.method){
			method.choice = FALSE
			
			if( sign.method=="limma"  & path.method=="MLP"  ){
				method.choice = TRUE
			}
			if(method.choice==TRUE){
				return(list(sign.method=sign.method,path.method=path.method))
			}	
			else{
				stop("Incorrect choice of method.")
			}
			
		}
		
		method.out = method.test(method[1],method[2])
		
		sign.method = method.out$sign.method
		path.method = method.out$path.method
		
		if(length(GeneInfo$ENTREZID)==1){
			GeneInfo$ENTREZID = colnames(GeneExpr)
		}
		
		# Determining the genesets if they were not given with the function input
		if((class(GENESET)=="geneSetMLP")[1] ){
			geneSet <- GENESET
		}
		else{
			geneSet <- MLP::getGeneSets(species = "Human",geneSetSource = geneSetSource,entrezIdentifiers = GeneInfo$ENTREZID)
		}
		
		if(is.null(topP)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}
		
		
		ResultMLP=list()
		
		for (k in 1:length(DataPrepared)){
			message(k)
					
			PathwaysResults=list()
	
			for (i in 1:length(DataPrepared$pvalsgenes)){
				message(paste(k,i,sep='.'))
				temp=list()
				temp[[1]]=DataPrepared[[k]]$Compounds[[i]] # the compounds
				temp[[2]]=DataPrepared[[k]]$Genes[[i]]		# the genes		
				pvalscluster=DataPrepared[[k]]$pvalsgenes[[i]]
				
				Entrezs=sapply(names(pvalscluster),function(x) return(GeneInfo$ENTREZID[which(GeneInfo$SYMBOL==x)]))
				
				if(path.method=="MLP"){
					## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
					
					names(pvalscluster) = Entrezs
					
					out.mlp <- MLP::MLP(
							geneSet = geneSet,
							geneStatistic = pvalscluster,
							minGenes = 5,
							maxGenes = 100,
							rowPermutations = TRUE,
							nPermutations = 100,
							smoothPValues = TRUE,
							probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9,addGeneSetDescription=TRUE)
					
					output = list()
					#output$gene.p.values = p.adjust(p.values,method="fdr")
					
					#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
					#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
					#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)
					
					#if(is.null(topP)){
					#		topP=length(ranked.genesets.table$p.values<=sign)
					#}
					
					if(is.null(topP)){
						topP=length(which(out.mlp$geneSetPValue<=sign))
					}
					
					#TopPaths=ranked.genesets.table[1:topP,]
					#AllPaths=ranked.genesets.table
					
					TopPaths=out.mlp[1:topP,]
					AllPaths=out.mlp
					
					output$TopPaths=TopPaths
					output$AllPaths=AllPaths
					#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]
					
					
					#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) 	)
					#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
					#output$nr.genesets = nr.genesets

					#output$object = out.mlp
					#output$method = "MLP"
					
					temp[[3]]=output				
				}
				names(temp)=c("Compounds","Genes","Pathways")
				PathwaysResults[[i]]=temp
				names(PathwaysResults)[i]=paste("Cluster",i,sep=" ")
				
			}
			
			ResultMLP[[k]]=PathwaysResults	
		}
		names(ResultMLP)=names
		for(i in 1:length(ResultMLP)){
			for(k in 1:length(ResultMLP[[i]])){
				if(is.null(ResultMLP[[i]][[k]])[1]){
					ResultMLP[[i]][[k]]=NA
					names(ResultMLP[[i]])[k]=paste("Cluster",k,sep=" ")
				}			
			}
			if(length(ResultMLP[[i]]) != maxclus){
				extra=maxclus-length(ResultMLP[[i]])
				for(j in 1:extra){
					ResultMLP[[i]][[length(ResultMLP[[i]])+j]]=NA
					names(ResultMLP[[i]])[length(ResultMLP[[i]])]=paste("Cluster",length(ResultMLP[[i]]),sep=" ")
				}
			}
		} 	
	}
	return(ResultMLP)	
}
