DiffGenesSelection=function(List,Selection,GeneExpr=NULL,nrclusters=NULL,method="limma",sign=0.05,topG=NULL,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	} 
	
	if(is.null(topG)){
		top1=FALSE
	}
	else{
		top1=TRUE
	}	
	
	if(class(Selection)=="character"){
		ResultLimma=list()
		Genes=list()
		temp=list()
		
		LeadCpds=Selection #names of the compounds
		OrderedCpds=colnames(GeneExpr)
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")
		
		label = rep(0,dim(GeneExpr)[2])
		label[which(colnames(GeneExpr)%in%Selection)] = 1
		label.factor = factor(label)
	
		
		if(class(GeneExpr)[1]=="ExpressionSet"){
			
			if (!requireNamespace("a4Base", quietly = TRUE)) {
				stop("a4Base needed for this function to work. Please install it.",
						call. = FALSE)
			}
			
			GeneExpr$LeadCmpds<-label.factor 
			DElead <- a4Base::limmaTwoLevels(GeneExpr,"LeadCpds")
			
			allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
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
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			
		}
		else{
			
			design = stats::model.matrix(~label.factor)
			fit = limma::lmFit(GeneExpr,design=design)
			fit = limma::eBayes(fit)
			
			allDE=limma::topTable(fit,coef=2,n=dim(GeneExpr)[1],adjust="fdr",sort.by="P")
			if(is.null(allDE$ID)){
				allDE$ID <- rownames(allDE)
			}
			else
			{
				allDE$ID=allDE$ID
			}
			if(top1==TRUE){
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			else if(top1==FALSE){
				topG=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}	
			
		}
		temp[[2]]=result
		
		names(temp)=c("Compounds","Genes")
		ResultLimma[[1]]=temp
		names(ResultLimma)="Selection"
		
	}
	else if(class(Selection)=="numeric" & !(is.null(List))){
	
		ListNew=list()
		element=0
		
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
		ResultLimma=list()
		for(k in 1:dim(Matrix)[1]){	
			cluster=Selection
			
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			
			Genes=list()
			temp=list()
			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the compounds
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")
			
			label = rep(0,dim(Matrix)[2])
			label[which(Matrix[k,]==cluster)] = 1
			label.factor = factor(label)
			
			GeneExpr.2=GeneExpr[,colnames(Matrix)]
			
			if(class(GeneExpr.2)[1]=="ExpressionSet"){
				
				if (!requireNamespace("a4Base", quietly = TRUE)) {
					stop("a4Base needed for this function to work. Please install it.",
							call. = FALSE)
				}
				
				GeneExpr.2$LeadCmpds<-label.factor 
				DElead <- a4Base::limmaTwoLevels(GeneExpr.2,"LeadCpds")
				
				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(top1==TRUE){
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}
				else if(top1==FALSE){
					topG=length(which(allDE$adj.P.Val<=sign))
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}
				
			}
			else{
				
				
				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(GeneExpr.2,design=design)
				fit = limma::eBayes(fit)
				
				allDE=limma::topTable(fit,coef=2,n=dim(GeneExpr)[1],adjust="fdr",sort.by="P")
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

					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}	
				
			}
			temp[[2]]=result
			
			names(temp)=c("Compounds","Genes")
			ResultLimma[[k]]=temp
			names(ResultLimma)[k]=paste(names[k],": Cluster", cluster, sep="")
		}		
	}

	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultLimma)
	
}
