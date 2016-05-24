DiffGenes=function(List,Selection=NULL,GeneExpr=NULL,nrclusters=NULL,method="limma",sign=0.05,topG=NULL,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	} 	
	if(!is.null(Selection)){
		ResultLimma=DiffGenesSelection(List,Selection,GeneExpr,nrclusters,method,sign,topG,fusionsLog,WeightClust,names)
	}
	else{

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
		
		if(is.null(names)){
			for(j in 1:length(List)){
				names[j]=paste("Method",j,sep=" ")	
			}
		}
		
		names(ListNew)=names
	
		if(is.null(topG)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}	
		
		
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		List=ListNew
		ResultLimma=list()
		maxclus=0
		for (k in 1:dim(MatrixClusters)[1]){
			clusters=MatrixClusters[k,]
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			Genes=list()
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			for (i in clust){
				
				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)] 
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds") #names of the compounds
				
				label = rep(0,length(names(clusters)))
				label[which(clusters==i)] = 1
				
				label.factor = factor(label)
				GeneExpr.2=GeneExpr[,names(clusters)]
				
				if(class(GeneExpr.2)[1]=="ExpressionSet"){
					
					if (!requireNamespace("a4Base", quietly = TRUE)) {
						stop("a4Base needed for this function to work. Please install it.",
								call. = FALSE)
					}
					
					
					GeneExpr.2$LeadCmpds<-label.factor	
					DElead <- a4Base::limmaTwoLevels(GeneExpr.2,"LeadCpds")
					
					allDE <-a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
					
					if(is.null(allDE$ID)){
						allDE$ID<- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					
					if(top1==TRUE){
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
						
					}
					else if(top1==FALSE){
						topG=length(which(allDE$adj.P.Val<=sign))
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
						
					}
					
				}
				else{

					design = stats::model.matrix(~label.factor)
					fit = limma::lmFit(GeneExpr.2,design=design)
					fit = limma::eBayes(fit)
					allDE=limma::topTable(fit,n=dim(GeneExpr)[1],coef=2,adjust="fdr",sort.by="P")
					
					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					
					if(top1==TRUE){
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
					}
					else if(top1==FALSE){
						topG=length(which(allDE$adj.P.Val<=sign))
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
					}
				}	
				
				temp[[2]]=result
				
				names(temp)=c("Compounds","Genes")
				
				Genes[[i]]=temp
				
				names(Genes)[i]=paste("Cluster",i,sep=" ")
			}
			ResultLimma[[k]]=Genes
			
		}
		names(ResultLimma)=names
		for(i in 1:length(ResultLimma)){
			for(k in 1:length(ResultLimma[[i]])){
				if(is.null(ResultLimma[[i]][[k]])[1]){
					ResultLimma[[i]][[k]]=NA
					names(ResultLimma[[i]])[k]=paste("Cluster",k,sep=" ")
				}			
			}
			if(length(ResultLimma[[i]]) != maxclus){
				extra=maxclus-length(ResultLimma[[i]])
				#temp=length(ResultLimma[[i]])
				for(j in 1:extra){
					ResultLimma[[i]][[length(ResultLimma[[i]])+1]]=NA
					names(ResultLimma[[i]])[length(ResultLimma[[i]])]=paste("Cluster",length(ResultLimma[[i]]),sep=" ")
				}
			}
		} 	

	}
	return(ResultLimma)
}
