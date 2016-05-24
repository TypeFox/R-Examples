ChooseCluster=function(Interactive=TRUE,LeadCpds=NULL,ClusterResult,ColorLab=NULL,BinData=NULL,ContData=NULL,Datanames=c("FP"),GeneExpr,topChar = 20, topG = 20,sign=0.05,nrclusters=NULL,cols=NULL,N=1){
	if(is.null(Datanames)){
		for(j in 1:(length(BinData)+length(ContData))){
			Datanames[j]=paste("Data",j,sep=" ")	
		}
	}
	OrInteractive=Interactive
	
	if(Interactive==TRUE){
		#windows()
		ClusterPlot(ClusterResult,ColorLab,nrclusters,cols)
		hc1<-stats::as.hclust(ClusterResult$Clust)
		ClusterSpecs<-list()
		ClusterSpecs=graphics::identify(hc1, N=N, MAXCLUSTER = nrow(BinData[[1]]), function(j) ChooseCluster(Interactive=FALSE,LeadCpds=rownames(BinData[[1]][j,]),ClusterResult,ColorLab=NULL,BinData,ContData,Datanames,GeneExpr,topChar,topG,sign,nrclusters,cols))		
		
		names(ClusterSpecs)<-sapply(seq(1,N),FUN=function(x) paste("Choice",x,sep=" "))
		
	}
	else{
		
		#BinData can be a list of different binary matrices :: different sources of information
		if(class(BinData)!="list"){
			stop("The binary data matrices must be put into a list")
		}
		for(i in 1:length(BinData)){
			BinData[[i]]=BinData[[i]]+0
			BinData[[i]]<-BinData[[i]][,which(colSums(BinData[[i]]) != 0 & colSums(BinData[[i]]) != nrow(BinData[[i]]))]
		}
		
		cpdSetG <-colnames(GeneExpr)
		
		
		DistM<-ClusterResult$DistM
		Clust<-ClusterResult$Clust
		if(is.null(Clust)){
			ClusterResult$Clust=ClusterResult
			Clust<-ClusterResult$Clust
		}
		
		hc <- stats::as.hclust(Clust)
		OrderedCpds <- hc$labels[hc$order]
		
		if(class(LeadCpds)=="character"){
			LeadCpds=list(LeadCpds)
		}
		
		Specs=list()
		for(i in 1:length(LeadCpds)){
			Compounds=list(LeadCpds[[i]],OrderedCpds)
			names(Compounds)=c("LeadCpds","OrderedCpds")
			cpdSet <- rownames(BinData[[1]])
			group <- factor(ifelse(cpdSet %in% LeadCpds[[i]], 1, 0)) #identify the group of interest

			#Determine characteristic features for the compounds: fishers exact test
			Characteristics=list()
			
			resultB=list()
			if(!is.null(BinData)){
			for(j in 1: length(BinData)){

				binMat=BinData[[j]]
				
				pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
				
				pFish <- sort(pFish)
				adjpFish<-stats::p.adjust(pFish, method = "fdr")
				
				AllFeat=data.frame(Names=names(pFish),P.Value=pFish,adj.P.Val=adjpFish)
				AllFeat$Names=as.character(AllFeat$Names)
				
				if(is.null(topChar)){
					topChar=length(which(pFish<sign))
				}
				
				TopFeat=AllFeat[0:topChar,]
				TopFeat$Names=as.character(TopFeat$Names)
				temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
				resultB[[j]]<-temp1
				names(resultB)[j]=Datanames[length(BinData)+j]
				
			}
			}
			resultC=list()
			if(!is.null(ContData)){
				for(j in 1:length(ContData)){
					contMat=ContData[[j]]
					
					group1=which(group==1)
					group2=which(group==0)
			
					pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)
					
					pTTest <- sort(pTTest)
					adjpTTest<-stats::p.adjust(pTTest, method = "fdr")
	
					AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topChar)){
						topC=length(which(pTTest<sign))
					}
					
					TopFeat=data.frame(Names=as.character(names(pTTest[0:topChar])),P.Value=pTTest[0:topChar],adj.P.Val=adjpTTest[0:topChar])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultC[[j]]<-temp1
					names(resultC)[j]=Datanames[j]
					
				}
			}
			Characteristics=c(resultB,resultC)
			names(Characteristics)=Datanames

			#Determine DE Genes with limma --> make difference between "regular" data matrix and "expression set"
			#GeneExpr.2=GeneExpr[,colnames(Matrix)]
			groupG <- factor(ifelse(cpdSetG %in% LeadCpds[[i]], 1, 0)) #identify the group of interest
			if(class(GeneExpr)[1]=="ExpressionSet"){
				GeneExpr$LeadCmpds<-groupG		
				
				
				if (!requireNamespace("a4Base", quietly = TRUE)) {
						stop("a4Base needed for this function to work. Please install it.",
								call. = FALSE)
				}

				DElead <- a4Base::limmaTwoLevels(GeneExpr,"LeadCpds")
				
				#allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), resort.by = "logFC",sort.by="p")
				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE <- allDE[0:topG, ]
				#TopAdjPval<-TopDE$adj.P.Val
				#TopRawPval<-TopDE$P.Value
				
				#RawpVal<-allDE$P.Value
				#AdjpVal <- allDE$adj.P.Val
				#genesEntrezId <- allDE$ENTREZID
				
				Genes<-list(TopDE,allDE)
				names(Genes)<-c("TopDE","AllDE")
				#Genes <- list(TopDE$SYMBOL,TopAdjPval,TopRawPval,genesEntrezId,RawpVal,AdjpVal)	
				#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
			}
			else{
				
				label.factor = factor(groupG)
				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(GeneExpr,design=design)
				fit = limma::eBayes(fit)
				
				#allDE = topTable(fit,coef=2,adjust="fdr",n=nrow(GeneExpr),resort.by = "logFC", sort.by="p")
				allDE = limma::topTable(fit,coef=2,adjust="fdr",n=nrow(GeneExpr), sort.by="p")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE=allDE[0:topG,]
				#TopAdjPval<-TopDE$adj.P.Val
				#TopRawPval<-TopDE$P.Value
				
				#RawpVal<-allDE$P.Value
				#AdjpVal <- allDE$adj.P.Val
				
				Genes<-list(TopDE,allDE)
				names(Genes)<-c("TopDE","AllDE")
				#Genes <- list(TopDE[,1],TopAdjPval,TopRawPval,allDE[,1],RawpVal,AdjpVal)	
				#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
				
			}
			
			out=list(Compounds,Characteristics,Genes)
			names(out)=c("Compounds","Characteristics","Genes")
			Specs[[i]]=out
			names(Specs)[i]=paste("Choice",i,sep=" ")
		}		
		if(OrInteractive==TRUE|length(Specs)==1){
			return(out)
		}
		else{
			return(Specs)
		}
	}
	class(ClusterSpecs)="ChosenClusters"
	return(ClusterSpecs)
}
