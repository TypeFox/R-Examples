PathwaysSelection<-function(List=NULL,Selection,GeneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),GeneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
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

		
	if(class(Selection)=="character"){
		ResultMLP=list()
		
		DataPrepared<-plyr::try_default(PreparePathway(List[[1]],GeneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(DataPrepared)){
			Temp=List[[1]]
			Compounds=list()
			Compounds$LeadCpds=Selection
			Compounds$OrderedCpds=colnames(GeneExpr)
			Temp[[length(Temp)+1]]=list(Compounds=Compounds)
			names(Temp)[length(Temp)]=paste("Cluster")
			
			DataPrepared<-PreparePathway(Temp,GeneExpr,topG,sign)
		}
		
		
		temp=list()
		temp[[1]]=DataPrepared$Compounds[[1]] #names of the compounds
		temp[[2]]=DataPrepared$Genes[[1]]
		
		
		
		
		if(path.method=="MLP"){
			## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
			p.values=DataPrepared$pvalsgenes[[1]]
			
			Entrezs=sapply(names(p.values),function(x) return(GeneInfo$ENTREZID[which(GeneInfo$SYMBOL==x)]))
			
			names(p.values) = Entrezs
			out.mlp <- MLP::MLP(
					geneSet = geneSet,
					geneStatistic = p.values,
					minGenes = 5,
					maxGenes = 100,
					rowPermutations = TRUE,
					nPermutations = 100,
					smoothPValues = TRUE,
					probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)
			
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
			
			#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
			#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
			#output$nr.genesets = nr.genesets
			
			#output$object = out.mlp
			#output$method = "MLP"
			
			temp[[3]]=output				
			
			
		}
		names(temp)=c("Compounds","Genes","Pathways")			
		ResultMLP[[1]]=temp
		names(ResultMLP)="Selection"
	}
	
		
	
	else if(class(Selection)=="numeric" & !(is.null(List))){
		check<-plyr::try_default(PreparePathway(List[[1]],GeneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(check)){
			
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
			Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
			List=ListNew	
			
			DataPrepared=list()
			for (k in 1:dim(Matrix)[1]){
							
				cluster=Selection
				
				check<-plyr::try_default(PreparePathway(List[[k]],GeneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(check)){
					Temp=List[[k]]
					Compounds=list()
					Compounds$LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)]
					Compounds$OrderedCpds=stats::as.hclust(List[[k]]$Clust$Clust)$labels[stats::as.hclust(List[[k]]$Clust$Clust)$order]
					Temp[[length(Temp)+1]]=list(Compounds=Compounds)
					names(Temp)[length(Temp)]=paste("Cluster")
					
					DataPrepared[[k]]<-PreparePathway(Temp,GeneExpr,topG,sign)
				}
			}
			names(DataPrepared)=names
		}
		
		else{
			for(k in 1:length(List)){
				DataPrepared[[k]]<-plyr::try_default(PreparePathway(List[[k]],GeneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(DataPrepared[[k]])){
					Temp=List[[k]]
					
				
					Compounds=list()
					Compounds$LeadCpds=List[[k]]$Compounds$LeadCpds
					Compounds$OrderedCpds=List[[k]]$Compounds$OrderedCpds
						
					Temp[[length(Temp)+1]]=list(Compounds=Compounds)
					names(Temp)[length(Temp)]=paste("Cluster")
						
					DataPrepared[[k]]<-PreparePathway(Temp,GeneExpr,topG,sign)
					
				}	
			}
			
			
		}
				
	ResultMLP=list()
	for (k in 1:length(DataPrepared)){
			message(k)
			cluster=Selection
			
			
			temp=list()
			temp[[1]]=DataPrepared[[k]]$Compounds[[1]] #names of the compounds
			temp[[2]]=DataPrepared[[k]]$Genes[[1]]
			
			
			if(path.method=="MLP"){
				## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
				p.values=DataPrepared[[k]]$pvalsgenes[[1]]
				Entrezs=sapply(names(p.values),function(x) return(GeneInfo$ENTREZID[which(GeneInfo$SYMBOL==x)]))
				
				names(p.values) = Entrezs
				
				out.mlp <- MLP::MLP(
						geneSet = geneSet,
						geneStatistic = p.values,
						minGenes = 5,
						maxGenes = 100,
						rowPermutations = TRUE,
						nPermutations = 100,
						smoothPValues = TRUE,
						probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)
				
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
				
				#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
				#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
				#output$nr.genesets = nr.genesets
				
				#output$object = out.mlp
				#output$method = "MLP"
				
				temp[[3]]=output				
			}
			names(temp)=c("Compounds","Genes","Pathways")			
			ResultMLP[[k]]=temp
			
		}
		names(ResultMLP)=names
				
	}
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	
	return(ResultMLP)	
}
